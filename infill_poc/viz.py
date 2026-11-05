import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

CMAP = plt.get_cmap("viridis")


def plot_voronoi_map(cells, xy, stiffness, infill, size_x, size_y, path, vmin=0, vmax=100):
    norm = Normalize(vmin, vmax)
    fig, ax = plt.subplots(figsize=(8, 8 * size_y / size_x + 0.6), dpi=200)
    for i, (cell, k, f) in enumerate(zip(cells, stiffness, infill)):
        cx, cy = cell.exterior.xy
        ax.fill(cx, cy, color=CMAP(norm(f)), edgecolor="black", linewidth=1.2)
        tx, ty = xy[i]
        ax.plot(tx, ty, "o", color="white", markeredgecolor="black", markersize=6)
        light = norm(f) < 0.55
        lx = float(np.clip(tx, 0.14 * size_x, 0.86 * size_x))
        above = ty < 0.25 * size_y
        ax.text(lx, ty + 1.4 if above else ty - 1.4, f"R{i}\nk={k:.2f} N/mm\n{f:.1f}%", ha="center",
                va="bottom" if above else "top", fontsize=8, color="white" if light else "black",
                fontweight="bold")
    ax.set_xlim(0, size_x)
    ax.set_ylim(0, size_y)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Voronoi regions from probe points, coloured by assigned infill")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="Assigned infill density (%)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_modifiers_3d(cells, infill, size_x, size_y, height, path, vmin=0, vmax=100):
    norm = Normalize(vmin, vmax)
    fig = plt.figure(figsize=(9, 6), dpi=200)
    ax = fig.add_subplot(projection="3d")
    for cell, f in zip(cells, infill):
        ring = np.asarray(cell.exterior.coords)[:-1]
        top = np.column_stack([ring, np.full(len(ring), height)])
        sides = [[(*ring[i], 0), (*ring[(i + 1) % len(ring)], 0), (*ring[(i + 1) % len(ring)], height), (*ring[i], height)]
                 for i in range(len(ring))]
        ax.add_collection3d(Poly3DCollection([top, *sides], facecolor=CMAP(norm(f)), edgecolor="black",
                                             linewidth=0.6))
    ax.set_xlim(0, size_x)
    ax.set_ylim(0, size_y)
    ax.set_zlim(0, height)
    ax.set_box_aspect((size_x, size_y, max(height, 0.35 * max(size_x, size_y))))
    ax.view_init(elev=30, azim=-55)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    ax.set_title(f"Modifier volumes extruded through full part height ({height:g} mm)")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.08, label="Infill density override (%)")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_verification(target, measured, tolerance, path):
    idx = np.arange(len(target))
    err = measured - target
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=200, gridspec_kw={"width_ratios": [1.6, 1]})
    a1.bar(idx - 0.2, target, 0.4, label="Target (from stiffness)", color="#9aa5b1")
    a1.bar(idx + 0.2, measured, 0.4, label="Measured from G-code", color="#1f6f8b")
    for i, e in zip(idx, err):
        a1.text(i, max(target[i], measured[i]) + 1.5, f"{e:+.1f}", ha="center", fontsize=8,
                color="black" if abs(e) <= tolerance else "#c0392b", fontweight="bold")
    a1.set_xticks(idx, [f"R{i}" for i in idx])
    a1.set_ylabel("Infill density (%)")
    a1.set_title("Per-region infill: target vs. measured (labels = error, pp)")
    a1.set_ylim(0, max(target.max(), measured.max()) + 20)
    a1.legend(loc="upper center", ncol=2)

    lim = (0, max(target.max(), measured.max()) + 5)
    a2.fill_between(lim, [lim[0] - tolerance, lim[1] - tolerance], [lim[0] + tolerance, lim[1] + tolerance],
                    color="#2ecc71", alpha=0.15, label=f"±{tolerance:g} pp tolerance")
    a2.plot(lim, lim, "k--", linewidth=1, label="y = x")
    a2.scatter(target, measured, color="#1f6f8b", zorder=3)
    for i in idx:
        a2.annotate(f"R{i}", (target[i], measured[i]), textcoords="offset points", xytext=(4, 4), fontsize=8)
    a2.set_xlim(lim)
    a2.set_ylim(lim)
    a2.set_aspect("equal")
    a2.set_xlabel("Target infill (%)")
    a2.set_ylabel("Measured infill (%)")
    a2.set_title("Agreement")
    a2.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_layer_toolpaths(cells, infill, infill_segs, offset_xy, size_x, size_y, path, vmin=0, vmax=100):
    """Actual sliced infill toolpaths for two consecutive mid-height layers over the Voronoi regions."""
    norm = Normalize(vmin, vmax)
    zs = np.unique(np.round(infill_segs[:, 4], 4))
    mid = len(zs) // 2
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.8), dpi=200)
    for ax, z in zip(axes, zs[mid:mid + 2]):
        for i, (cell, f) in enumerate(zip(cells, infill)):
            cx, cy = cell.exterior.xy
            ax.fill(cx, cy, color=CMAP(norm(f)), alpha=0.35, edgecolor="none")
            c = cell.representative_point()
            ax.text(c.x, c.y, f"R{i}\n{f:.0f}%", ha="center", va="center", fontsize=8, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
        layer = infill_segs[np.isclose(infill_segs[:, 4], z)]
        for s in layer:
            ax.plot([s[0] - offset_xy[0], s[2] - offset_xy[0]], [s[1] - offset_xy[1], s[3] - offset_xy[1]],
                    color="black", linewidth=0.6)
        for cell in cells:
            ax.plot(*cell.exterior.xy, color="#c0392b", linewidth=1.4)
        ax.set_xlim(0, size_x)
        ax.set_ylim(0, size_y)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"Sliced infill toolpaths, layer z = {z:.2f} mm")
    fig.suptitle("PrusaSlicer output: infill line spacing changes at each Voronoi boundary (red)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
