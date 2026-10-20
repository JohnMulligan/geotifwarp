# geotifwarp

Warps a scanned map onto modern coordinates using the control points you placed in the QGIS Georeferencer, and writes the result as a standard XYZ tile pyramid (`{z}/{x}/{y}.png`, Web Mercator, transparent outside the map).

It uses a piecewise-affine (triangulation) rubber-sheet: every enabled point lands exactly where you put it, and each point only affects the triangles that touch it, so adding a point can't disturb distant parts of the map.

## Setup

```sh
cd ~/code/geotifwarp
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 1. Export the points from QGIS

In the Georeferencer, open the scan, place your points, then **File → Save GCP Points as…** and save `map.tif.points`.

- Unticked points are ignored.
- Make sure the CRS QGIS asks for when saving is the CRS you actually picked the points in. It's written to the file's `#CRS:` line, and the script reads it from there.
- You don't need to run the Georeferencer's own transform.

## 2. Build the tiles

```sh
python geotifwarp.py map.tif map.tif.points map_tiles/
```

This produces:

```
map_tiles/
  metadata.json        bounds, zoom range, tile counts
  12/1369/1897.png
  ...
  17/43835/60713.png
```

Zoom levels are chosen automatically: the top level is the first one at or finer than the scan's resolution, and the bottom level is where the whole map fits in about one tile. Tiles with nothing on them aren't written. The output folder has to be new or empty, so a re-run never mixes old and new tiles.

| Option | Use |
|---|---|
| `--min-zoom N` / `--max-zoom N` | Override the zoom range. Each extra zoom level at the top roughly quadruples the tile count. |
| `--resampling nearest\|bilinear\|cubic` | Default is `bilinear`. Paletted scans always use `nearest`. |
| `--workers N` | Parallel processes. Defaults to all cores. |
| `--points-crs EPSG:xxxx` | Only needed if the `.points` file has no `#CRS:` line (older QGIS). |

Input must be an 8-bit scan: greyscale, RGB, RGBA or paletted.

## 3. View it in QGIS

**Layer → Add Layer → Add XYZ Layer… → New**. Set the URL to

```
file:///Users/johnmulligan/path/to/map_tiles/{z}/{x}/{y}.png
```

and set the max zoom to the `maxzoom` value in `metadata.json`. The same folder can be served as-is to Leaflet, OpenLayers or MapLibre.

## Placing points for this method

- **Only the area inside your points is output.** Anything outside the outer ring of points is transparent, so put points at the corners and along the edges of the area you want to keep.
- **A "mirrored triangles" warning** means two or more points are in a different relative order on the scan than on the map, which makes the image fold over itself. The script lists the map coordinates involved. Usually one of them is misplaced or misidentified.
- **Creases at triangle edges** are expected where neighbouring points pull in very different directions. Adding a point in the middle of the creased area smooths it out.
- Scans you load in QGIS without georeferencing and scans that already had rough georeferencing both work. The script detects which and reads the source coordinates accordingly.
