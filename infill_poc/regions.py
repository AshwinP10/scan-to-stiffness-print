import numpy as np
import shapely
import trimesh
from shapely.geometry import MultiPoint, Point, box


def make_box_model(size_x, size_y, height):
    """Axis-aligned test part spanning [0, sx] x [0, sy] x [0, h] in model coordinates."""
    mesh = trimesh.creation.box(extents=(size_x, size_y, height))
    mesh.apply_translation((size_x / 2, size_y / 2, height / 2))
    return mesh


def sample_probe_points(n, size_x, size_y, k_range, min_spacing, rng, edge_margin=2.0, max_tries=20000):
    """Random top-face probe locations with a minimum spacing, plus a random stiffness per point.

    The spacing floor keeps every Voronoi cell wide enough to hold several infill lines;
    at 10% density the line pitch is ~4 mm, so tiny cells would be unmeasurable.
    """
    pts = []
    for _ in range(max_tries):
        p = rng.uniform((edge_margin, edge_margin), (size_x - edge_margin, size_y - edge_margin))
        if all(np.hypot(*(p - q)) >= min_spacing for q in pts):
            pts.append(p)
            if len(pts) == n:
                break
    if len(pts) < n:
        raise ValueError(f"Could only place {len(pts)}/{n} points with min_spacing={min_spacing} mm; "
                         "lower --n-points or --min-spacing.")
    xy = np.array(pts)
    k = rng.uniform(k_range[0], k_range[1], n)
    return xy, k


def voronoi_cells(xy, bounds):
    """One clipped, convex Voronoi cell per input point, in input order."""
    cells = shapely.voronoi_polygons(MultiPoint(xy), extend_to=bounds, ordered=True)
    out = []
    for p, cell in zip(xy, cells.geoms):
        clipped = cell.intersection(bounds).simplify(1e-9)
        if not clipped.contains(Point(p)):
            raise RuntimeError("Voronoi cell order does not match input points")
        out.append(clipped)
    return out


def prism_mesh(poly, z0, z1):
    """Watertight extrusion of a convex polygon between z0 and z1, outward-facing normals."""
    if abs(poly.convex_hull.area - poly.area) > 1e-6 * max(poly.area, 1.0):
        raise ValueError("prism_mesh expects a convex polygon")
    ring = shapely.geometry.polygon.orient(poly, sign=1.0).exterior.coords[:-1]
    ring = np.asarray(ring)
    n = len(ring)
    verts = np.vstack([np.column_stack([ring, np.full(n, z0)]),
                       np.column_stack([ring, np.full(n, z1)])])
    faces = []
    for i in range(1, n - 1):
        faces.append([0, i + 1, i])
        faces.append([n, n + i, n + i + 1])
    for i in range(n):
        j = (i + 1) % n
        faces.append([i, j, n + j])
        faces.append([i, n + j, n + i])
    mesh = trimesh.Trimesh(verts, np.array(faces), process=False)
    if not mesh.is_watertight or mesh.volume <= 0:
        raise RuntimeError("Generated modifier prism is not a valid closed solid")
    return mesh


def build_regions(xy, size_x, size_y, height, wall_overhang=0.5, seam_overlap=0.05):
    """Returns (measurement_cells, modifier_meshes).

    measurement_cells are the exact Voronoi cells clipped to the part footprint (used for G-code
    binning). Modifiers overhang the outer walls and overlap their neighbours very slightly so no
    sliver along a coincident seam falls back to the base infill.
    """
    footprint = box(0, 0, size_x, size_y)
    extended = box(-wall_overhang, -wall_overhang, size_x + wall_overhang, size_y + wall_overhang)
    ext_cells = voronoi_cells(xy, extended)
    measurement_cells = [c.intersection(footprint) for c in ext_cells]
    modifiers = []
    for c in ext_cells:
        grown = c.buffer(seam_overlap, join_style="mitre").intersection(extended).simplify(1e-9)
        modifiers.append(prism_mesh(grown, 0.0, height))
    return measurement_cells, modifiers
