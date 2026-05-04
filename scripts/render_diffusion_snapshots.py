#!/usr/bin/env python3
"""Render heterogeneous-D diffusion demo snapshots to PNG.

Reads legacy ASCII VTK STRUCTURED_POINTS files written by
``examples/heterogeneous_diffusion.cpp`` (via
``examples/common/vtk_writer.hpp``) and emits one PNG per input file
showing a 2D z-midplane slice of psi.

Usage::

    python3 scripts/render_diffusion_snapshots.py \
        --in-dir /tmp/hd_figs \
        --out-dir docs/user-guide/figures/heterogeneous-diffusion \
        --steps 200 1000 4000

Each output PNG is named ``early.png``, ``mid.png``, ``late.png`` by
position in ``--steps`` (or ``snapshot_<step>.png`` if more than three
steps are passed). The demo's diffusivity varies smoothly in `x`, so
the Gaussian spreads anisotropically: faster in high-D regions
(near `x = 0` and `x = L`), slower in low-D regions (near
`x = L/2`). Mirrors the Phase 12 ``render_ch_snapshots.py`` pattern
with an adjusted sequential colormap (positive Gaussian, not a
±1-bounded order parameter).

Requires ``numpy`` and ``matplotlib``.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np


_NAME_BY_INDEX = ("early", "mid", "late")


def parse_vtk(path: Path) -> tuple[tuple[int, int, int], np.ndarray]:
    """Parse a legacy ASCII VTK STRUCTURED_POINTS file.

    Returns ``((Nx, Ny, Nz), data)`` with ``data`` shape ``(Nz, Ny, Nx)``
    so that ``data[k, :, :]`` is a 2D slice at z-index ``k``.
    """
    with path.open() as fh:
        lines = fh.read().splitlines()

    dims_line = next(line for line in lines if line.startswith("DIMENSIONS "))
    Nx, Ny, Nz = (int(t) for t in dims_line.split()[1:4])

    data_start = lines.index("LOOKUP_TABLE default") + 1
    flat = np.fromstring("\n".join(lines[data_start:]), sep="\n", dtype=np.float64)
    if flat.size != Nx * Ny * Nz:
        raise ValueError(
            f"{path}: expected {Nx * Ny * Nz} values, got {flat.size}"
        )
    return (Nx, Ny, Nz), flat.reshape((Nz, Ny, Nx))


def render_slice(field_3d: np.ndarray,
                 title: str,
                 out_path: Path,
                 vmax: float | None) -> None:
    """Render the z-midplane slice as a PNG with a sequential colormap.

    The diffusion demo's psi is non-negative, so a sequential colormap
    (``viridis``) reads better than the diverging ``RdBu_r`` used by
    the CH demo.
    """
    Nz = field_3d.shape[0]
    slice_2d = field_3d[Nz // 2]
    if vmax is None:
        vmax = max(1e-6, float(slice_2d.max()))

    fig, ax = plt.subplots(figsize=(4, 4))
    img = ax.imshow(
        slice_2d,
        origin="lower",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    fig.colorbar(img, ax=ax, fraction=0.046, pad=0.04, label=r"$\psi$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def find_snapshot(in_dir: Path, step: int) -> Path:
    pattern = re.compile(rf"psi_0*{step}\.vtk$")
    for path in sorted(in_dir.iterdir()):
        if pattern.search(path.name):
            return path
    raise FileNotFoundError(f"no psi_*{step:06d}.vtk under {in_dir}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--steps", type=int, nargs="+", required=True)
    parser.add_argument("--shared-vmax", action="store_true",
                        help="Use the same colormap range across all output "
                             "PNGs (= max value found in the first / earliest "
                             "snapshot). Default: per-snapshot autoscale.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    snapshots = []
    for step in args.steps:
        snap = find_snapshot(args.in_dir, step)
        dims, data = parse_vtk(snap)
        snapshots.append((step, snap, dims, data))

    shared_vmax: float | None = None
    if args.shared_vmax:
        Nz = snapshots[0][3].shape[0]
        shared_vmax = float(snapshots[0][3][Nz // 2].max())

    for idx, (step, snap, dims, data) in enumerate(snapshots):
        if idx < len(_NAME_BY_INDEX) and len(snapshots) <= len(_NAME_BY_INDEX):
            out_name = f"{_NAME_BY_INDEX[idx]}.png"
        else:
            out_name = f"snapshot_{step:06d}.png"
        out_path = args.out_dir / out_name
        title = f"step={step}"
        render_slice(data, title, out_path, shared_vmax)
        Nx, Ny, Nz = dims
        print(f"{snap.name} ({Nx}x{Ny}x{Nz}) -> {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
