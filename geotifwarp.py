#!/usr/bin/env python3
import argparse
import atexit
import json
import math
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import rasterio
import rasterio.errors
from PIL import Image
from rasterio.crs import CRS
from rasterio.warp import transform as reproject_points
from rasterio.windows import Window
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import map_coordinates

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)

RESAMPLING_ORDER = {"nearest": 0, "bilinear": 1, "cubic": 3}
WEB_MERCATOR = CRS.from_epsg(3857)
WORLD = 2 * math.pi * 6378137.0
HALF = WORLD / 2
TILE = 256

_state = {}


def read_points(path):
    crs = None
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            if line.startswith("#CRS:"):
                wkt = line[len("#CRS:"):].strip()
                if wkt:
                    crs = CRS.from_user_input(wkt)
                break
            if not line.startswith("#"):
                break
    df = pd.read_csv(path, comment="#", encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    if "enable" in df.columns:
        df = df[df["enable"].astype(int) == 1]
    return df, crs


def source_pixels(df, src):
    sx = df["sourceX"].to_numpy(float)
    sy = df["sourceY"].to_numpy(float)
    if src.transform.is_identity:
        return sx, -sy
    cols, rows = ~src.transform * (sx, sy)
    return np.asarray(cols), np.asarray(rows)


def signed_areas(simplices, p):
    p0, p1, p2 = p[simplices[:, 0]], p[simplices[:, 1]], p[simplices[:, 2]]
    return (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])


def tile_res(z):
    return WORLD / TILE / 2 ** z


def tile_range(z, minx, miny, maxx, maxy):
    n = 2 ** z
    size = WORLD / n
    x0 = max(int(math.floor((minx + HALF) / size)), 0)
    x1 = min(int(math.floor((maxx + HALF) / size)), n - 1)
    y0 = max(int(math.floor((HALF - maxy) / size)), 0)
    y1 = min(int(math.floor((HALF - miny) / size)), n - 1)
    return x0, x1, y0, y1


def tile_path(out, z, x, y):
    return os.path.join(out, str(z), str(x), f"{y}.png")


def save_tile(out, z, x, y, rgba):
    path = tile_path(out, z, x, y)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    Image.fromarray(rgba, "RGBA").save(path, optimize=False)


def init_worker(input_path, dst_xy, src_cr, order, out):
    src = rasterio.open(input_path)
    atexit.register(src.close)
    colormap = None
    try:
        colormap = src.colormap(1)
    except ValueError:
        pass
    lut = None
    if colormap:
        lut = np.zeros((256, 3), np.uint8)
        for k, v in colormap.items():
            if 0 <= k < 256:
                lut[k] = v[:3]
    _state.update(src=src, interp=LinearNDInterpolator(dst_xy, src_cr), order=0 if lut is not None else order, lut=lut, out=out)


def render_block(job):
    z, tx0, ty0, ntx, nty, tx_max, ty_max = job
    src, interp, order, lut, out = _state["src"], _state["interp"], _state["order"], _state["lut"], _state["out"]
    res = tile_res(z)
    w, h = ntx * TILE, nty * TILE
    left = tx0 * TILE * res - HALF
    top = HALF - ty0 * TILE * res
    gx, gy = np.meshgrid(left + (np.arange(w) + 0.5) * res, top - (np.arange(h) + 0.5) * res)
    cr = interp(gx.ravel(), gy.ravel())
    c = cr[:, 0] - 0.5
    r = cr[:, 1] - 0.5
    valid = ~np.isnan(c) & (c >= -0.5) & (c <= src.width - 0.5) & (r >= -0.5) & (r <= src.height - 0.5)
    if not valid.any():
        return []

    bands = [1] if lut is not None else ([1, 2, 3] if src.count >= 3 else [1])
    cv, rv = c[valid], r[valid]
    pad = 2 + order
    c0 = max(int(math.floor(cv.min())) - pad, 0)
    r0 = max(int(math.floor(rv.min())) - pad, 0)
    c1 = min(int(math.ceil(cv.max())) + pad + 1, src.width)
    r1 = min(int(math.ceil(rv.max())) + pad + 1, src.height)
    swin = Window(c0, r0, c1 - c0, r1 - r0)
    block = src.read(bands, window=swin)
    smask = src.dataset_mask(window=swin)
    coords = np.vstack([rv - r0, cv - c0])

    sampled = []
    for b in range(len(bands)):
        vals = map_coordinates(block[b].astype(np.float32), coords, order=order, mode="nearest", prefilter=order > 1)
        sampled.append(np.clip(np.rint(vals), 0, 255).astype(np.uint8))
    if lut is not None:
        rgb = lut[sampled[0]]
    elif len(sampled) == 1:
        rgb = np.repeat(sampled[0][:, None], 3, axis=1)
    else:
        rgb = np.stack(sampled, axis=1)

    rgba = np.zeros((h * w, 4), np.uint8)
    rgba[valid, :3] = rgb
    rgba[valid, 3] = np.where(map_coordinates(smask, coords, order=0, mode="nearest") > 0, 255, 0)
    rgba = rgba.reshape(h, w, 4)

    written = []
    for j in range(nty):
        for i in range(ntx):
            x, y = tx0 + i, ty0 + j
            if x > tx_max or y > ty_max:
                continue
            tile = rgba[j * TILE:(j + 1) * TILE, i * TILE:(i + 1) * TILE]
            if tile[..., 3].any():
                save_tile(out, z, x, y, np.ascontiguousarray(tile))
                written.append((x, y))
    return written


def render_parent(job):
    z, x, y, out = job
    canvas = Image.new("RGBA", (TILE * 2, TILE * 2), (0, 0, 0, 0))
    for dx in (0, 1):
        for dy in (0, 1):
            p = tile_path(out, z + 1, 2 * x + dx, 2 * y + dy)
            if os.path.exists(p):
                with Image.open(p) as child:
                    canvas.paste(child.convert("RGBA"), (dx * TILE, dy * TILE))
    small = canvas.resize((TILE, TILE), Image.Resampling.BOX)
    save_tile(out, z, x, y, np.asarray(small))
    return x, y


def main():
    ap = argparse.ArgumentParser(description="Piecewise-affine (triangulation) warp of a scanned map into an XYZ tile pyramid, using QGIS .points GCPs.")
    ap.add_argument("input", help="the scan you georeferenced (any raster GDAL reads)")
    ap.add_argument("points", help="QGIS Georeferencer .points file")
    ap.add_argument("output", help="output folder for {z}/{x}/{y}.png tiles (must be new or empty)")
    ap.add_argument("--points-crs", help="CRS of mapX/mapY when the .points file has no #CRS line, e.g. EPSG:32621")
    ap.add_argument("--min-zoom", type=int, help="lowest zoom to build (default: the map fits in about one tile)")
    ap.add_argument("--max-zoom", type=int, help="highest zoom to build (default: first zoom at or finer than the scan's resolution)")
    ap.add_argument("--resampling", choices=RESAMPLING_ORDER, default="bilinear")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()

    if os.path.isdir(args.output) and os.listdir(args.output):
        sys.exit(f"{args.output} is not empty; choose a new folder so old tiles don't mix with new ones.")

    df, pts_crs = read_points(args.points)
    if args.points_crs:
        pts_crs = CRS.from_user_input(args.points_crs)
    if pts_crs is None:
        sys.exit("No #CRS line in the .points file; pass --points-crs.")
    if len(df) < 3:
        sys.exit(f"Need at least 3 enabled points, found {len(df)}.")

    with rasterio.open(args.input) as src:
        if any(dt != "uint8" for dt in src.dtypes):
            sys.exit(f"Only 8-bit scans are supported (this one is {src.dtypes[0]}).")
        cols, rows = source_pixels(df, src)

    mx, my = reproject_points(pts_crs, WEB_MERCATOR, df["mapX"].to_numpy(float), df["mapY"].to_numpy(float))
    dst_xy = np.column_stack([mx, my])
    src_cr = np.column_stack([cols, rows])
    _, keep = np.unique(dst_xy.round(3), axis=0, return_index=True)
    if len(keep) < len(dst_xy):
        print(f"warning: dropped {len(dst_xy) - len(keep)} duplicate map coordinates", file=sys.stderr)
        keep.sort()
        dst_xy, src_cr = dst_xy[keep], src_cr[keep]

    tri = LinearNDInterpolator(dst_xy, src_cr).tri
    folds = tri.simplices[np.sign(signed_areas(tri.simplices, dst_xy)) != -np.sign(signed_areas(tri.simplices, src_cr))]
    if len(folds):
        orig = np.column_stack([df["mapX"].to_numpy(float), df["mapY"].to_numpy(float)])
        if len(keep) < len(orig):
            orig = orig[keep]
        print(f"warning: {len(folds)} triangle(s) are mirrored between map and scan, so those areas will fold. "
              "Check the order of these points (mapX, mapY):", file=sys.stderr)
        for i in sorted(set(folds.ravel())):
            print(f"  {orig[i, 0]:.3f}, {orig[i, 1]:.3f}", file=sys.stderr)

    native = math.sqrt(np.abs(signed_areas(tri.simplices, dst_xy)).sum() / np.abs(signed_areas(tri.simplices, src_cr)).sum())
    minx, miny = dst_xy.min(axis=0)
    maxx, maxy = dst_xy.max(axis=0)
    max_zoom = args.max_zoom if args.max_zoom is not None else min(max(math.ceil(math.log2(WORLD / TILE / native)), 0), 24)
    min_zoom = args.min_zoom if args.min_zoom is not None else min(max(int(math.floor(math.log2(WORLD / max(maxx - minx, maxy - miny)))), 0), max_zoom)
    if min_zoom > max_zoom:
        sys.exit("--min-zoom is greater than --max-zoom.")

    print(f"{len(dst_xy)} points, {len(tri.simplices)} triangles; scan is ~{native:.3g} m/px (Web Mercator), "
          f"building zoom {min_zoom}-{max_zoom}", file=sys.stderr)

    os.makedirs(args.output, exist_ok=True)
    x0, x1, y0, y1 = tile_range(max_zoom, minx, miny, maxx, maxy)
    step = 4
    jobs = [(max_zoom, tx, ty, step, step, x1, y1) for ty in range(y0, y1 + 1, step) for tx in range(x0, x1 + 1, step)]

    init = (args.input, dst_xy, src_cr, RESAMPLING_ORDER[args.resampling], args.output)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=init) as pool:
        present = set()
        for n, written in enumerate(pool.map(render_block, jobs), 1):
            present.update(written)
            print(f"\r  zoom {max_zoom}: block {n}/{len(jobs)}", end="", file=sys.stderr)
        print(f"\r  zoom {max_zoom}: {len(present)} tiles" + " " * 20, file=sys.stderr)
        counts = {max_zoom: len(present)}

        for z in range(max_zoom - 1, min_zoom - 1, -1):
            parents = sorted({(x // 2, y // 2) for x, y in present})
            present = set(pool.map(render_parent, [(z, x, y, args.output) for x, y in parents], chunksize=16))
            counts[z] = len(present)
            print(f"  zoom {z}: {len(present)} tiles", file=sys.stderr)

    (lon0, lon1), (lat0, lat1) = reproject_points(WEB_MERCATOR, CRS.from_epsg(4326), [minx, maxx], [miny, maxy])
    meta = {
        "name": os.path.splitext(os.path.basename(args.input))[0],
        "format": "png", "scheme": "xyz", "minzoom": min_zoom, "maxzoom": max_zoom,
        "bounds": [lon0, lat0, lon1, lat1], "center": [(lon0 + lon1) / 2, (lat0 + lat1) / 2, min_zoom],
        "tiles": counts,
    }
    with open(os.path.join(args.output, "metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"wrote {sum(counts.values())} tiles to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
