"""Measure the infill density PrusaSlicer actually produced inside each Voronoi region."""
import math
import re

import numpy as np
import shapely
from shapely.geometry import LineString, box

_WORD = re.compile(r"([A-Z])(-?\d*\.?\d+)")
INFILL = "Internal infill"
EXT_PERIMETER = "External perimeter"


def parse_extrusions(gcode_path, types=(INFILL, EXT_PERIMETER)):
    """Return {type: array[(x0, y0, x1, y1, z, e)]} for every extruding XY move of the given types.
    Handles both relative (M83) and absolute (M82 + G92) extrusion."""
    x = y = z = 0.0
    e_last = 0.0
    relative_e = False
    cur_type = None
    out = {t: [] for t in types}
    with open(gcode_path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            if raw.startswith(";TYPE:"):
                cur_type = raw[6:].strip()
                continue
            code = raw.split(";", 1)[0].strip().upper()
            if not code:
                continue
            w = dict(_WORD.findall(code))
            g, m = w.get("G"), w.get("M")
            if m == "83":
                relative_e = True
            elif m == "82":
                relative_e = False
            elif g == "92" and "E" in w:
                e_last = float(w["E"])
            elif g in ("0", "1"):
                nx, ny, nz = float(w.get("X", x)), float(w.get("Y", y)), float(w.get("Z", z))
                e = 0.0
                if "E" in w:
                    ev = float(w["E"])
                    e = ev if relative_e else ev - e_last
                    e_last = e_last if relative_e else ev
                if e > 0 and (nx != x or ny != y) and cur_type in out:
                    out[cur_type].append((x, y, nx, ny, nz, e))
                x, y, z = nx, ny, nz
    return {t: np.array(v, dtype=float).reshape(-1, 6) for t, v in out.items()}


def measured_offset(ext_perimeter_segs, size_x, size_y):
    """Bed-to-model XY offset, from the centre of the external perimeter's bounding box."""
    xs = np.concatenate([ext_perimeter_segs[:, 0], ext_perimeter_segs[:, 2]])
    ys = np.concatenate([ext_perimeter_segs[:, 1], ext_perimeter_segs[:, 3]])
    return np.array([(xs.min() + xs.max()) / 2 - size_x / 2, (ys.min() + ys.max()) / 2 - size_y / 2])


def region_densities(infill_segs, cells, offset_xy, layer_height, filament_diameter):
    """Per-region infill density (%) = deposited infill area / (fillable area x infill layers).

    A rectilinear line of length L at flow spacing s deposits L*h*s of plastic, and PrusaSlicer
    sets the line pitch to s / density, so (extruded volume / h) / fillable area recovers density.
    """
    seg = infill_segs.copy()
    seg[:, [0, 2]] -= offset_xy[0]
    seg[:, [1, 3]] -= offset_xy[1]
    n_layers = len(np.unique(np.round(seg[:, 4], 4)))

    # Infill centrelines are clipped to PrusaSlicer's fill polygon, so their extent is that polygon.
    xs = np.concatenate([seg[:, 0], seg[:, 2]])
    ys = np.concatenate([seg[:, 1], seg[:, 3]])
    fill_rect = box(xs.min(), ys.min(), xs.max(), ys.max())

    area_per_e = math.pi * (filament_diameter / 2) ** 2 / layer_height
    deposited = np.zeros(len(cells))
    in_a = np.stack([shapely.intersects_xy(c, seg[:, 0], seg[:, 1]) for c in cells], axis=1)
    in_b = np.stack([shapely.intersects_xy(c, seg[:, 2], seg[:, 3]) for c in cells], axis=1)
    clean = (in_a.sum(1) == 1) & (in_b.sum(1) == 1) & (in_a & in_b).any(1)
    owner = np.argmax(in_a, axis=1)
    np.add.at(deposited, owner[clean], seg[clean, 5] * area_per_e)

    n_split = 0
    for s in seg[~clean]:
        line = LineString([(s[0], s[1]), (s[2], s[3])])
        for i, c in enumerate(cells):
            frac = line.intersection(c).length / line.length
            if frac > 0:
                deposited[i] += s[5] * area_per_e * frac
        n_split += 1

    fill_areas = np.array([c.intersection(fill_rect).area for c in cells])
    density = 100.0 * deposited / (fill_areas * n_layers)
    return density, {"n_infill_layers": n_layers, "n_segments": len(seg), "n_split_segments": n_split,
                     "fill_rect": fill_rect.bounds, "fill_areas": fill_areas}
