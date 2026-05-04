#!/usr/bin/env python3
"""Render Cahn--Hilliard demo snapshots to PNG.

Reads legacy ASCII VTK STRUCTURED_POINTS files written by
``examples/cahn_hilliard.cpp`` (via ``examples/common/vtk_writer.hpp``)
and emits one PNG per input file showing a 2D z-midplane slice of the
order parameter psi.

Usage::

    python3 scripts/render_ch_snapshots.py \
        --in-dir /tmp/ch_figs \
        --out-dir docs/user-guide/figures/cahn-hilliard \
        --steps 1500 4000 8000

Each output PNG is named ``early.png``, ``mid.png``, ``late.png`` by
position in ``--steps`` (or ``snapshot_<step>.png`` if more than three
steps are passed). Requires ``numpy`` and ``matplotlib``.
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

    Returns a ``((Nx, Ny, Nz), data)`` tuple, where ``data`` has shape
    ``(Nz, Ny, Nx)`` (z-major) so a slice ``data[k, :, :]`` is a 2D
    image at z-index ``k``.
    """
    with path.open() as fh:
        lines = fh.read().splitlines()

    # Header lines are positional in the STRUCTURED_POINTS legacy format
    # written by writeVtkStructuredPoints; the data values follow once
    # we hit "LOOKUP_TABLE default".
    dims_line = next(line for line in lines if line.startswith("DIMENSIONS "))
    Nx, Ny, Nz = (int(t) for t in dims_line.split()[1:4])

    # Values are emitted one per line in i-fastest order (matches the
    # Field<double> storage layout `linear_index = i + Nx * (j + Ny * k)`).
    data_start = lines.index("LOOKUP_TABLE default") + 1
    flat = np.fromstring("\n".join(lines[data_start:]), sep="\n", dtype=np.float64)
    if flat.size != Nx * Ny * Nz:
        raise ValueError(
            f"{path}: expected {Nx * Ny * Nz} values, got {flat.size}"
        )
    # Reshape: i varies fastest, then j, then k.
    return (Nx, Ny, Nz), flat.reshape((Nz, Ny, Nx))


def render_slice(field_3d: np.ndarray, title: str, out_path: Path) -> None:
    """Render the z-midplane slice as a PNG with a diverging colormap."""
    Nz = field_3d.shape[0]
    slice_2d = field_3d[Nz // 2]
    vmax = max(0.05, float(np.abs(slice_2d).max()))

    fig, ax = plt.subplots(figsize=(4, 4))
    img = ax.imshow(
        slice_2d,
        origin="lower",
        cmap="RdBu_r",
        vmin=-vmax,
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
    parser.add_argument("--in-dir", type=Path, required=True,
                        help="Directory containing psi_NNNNNN.vtk snapshots.")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Directory to write the PNGs to.")
    parser.add_argument("--steps", type=int, nargs="+", required=True,
                        help="Step numbers to render. First three are named "
                             "early.png, mid.png, late.png; further steps use "
                             "snapshot_<step>.png.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    for idx, step in enumerate(args.steps):
        snap = find_snapshot(args.in_dir, step)
        (Nx, Ny, Nz), data = parse_vtk(snap)
        if idx < len(_NAME_BY_INDEX) and len(args.steps) <= len(_NAME_BY_INDEX):
            out_name = f"{_NAME_BY_INDEX[idx]}.png"
        else:
            out_name = f"snapshot_{step:06d}.png"
        out_path = args.out_dir / out_name
        title = f"step={step}"
        render_slice(data, title, out_path)
        print(f"{snap.name} ({Nx}x{Ny}x{Nz}) -> {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
