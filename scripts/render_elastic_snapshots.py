#!/usr/bin/env python3
"""Render elastic-energy demo snapshots to PNG.

Reads legacy ASCII VTK STRUCTURED_POINTS files written by
``examples/elastic_energy.cpp`` (via ``examples/common/vtk_writer.hpp``)
and emits one PNG per input file showing a 2D y-midplane slice of the
axial-strain field epsilon_axis,axis.

Usage::

    # produce three runs with different alpha amplitudes
    for alpha_label in 0.005 0.01 0.02; do
        ./build/clang-debug/examples/elastic_energy --n-x 64 --alpha $alpha_label \
            --out-dir /tmp/elastic_figs/$alpha_label
        cp /tmp/elastic_figs/$alpha_label/elastic_energy.vtk \
            /tmp/elastic_figs/elastic_alpha_${alpha_label}.vtk
    done
    python3 scripts/render_elastic_snapshots.py \\
        --in-dir /tmp/elastic_figs \\
        --out-dir docs/user-guide/figures/elastic-energy \\
        --files elastic_alpha_0.005.vtk elastic_alpha_0.01.vtk elastic_alpha_0.02.vtk

The output PNGs are named ``small.png``, ``medium.png``, ``large.png``
by position in ``--files`` (or ``snapshot_<basename>.png`` if more than
three are passed). The strain profile is `alpha * cos(k0 x_axis)` so
each PNG shows a single sin-wave-derivative pattern across the
midplane.

Requires ``numpy`` and ``matplotlib``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless rendering
import matplotlib.pyplot as plt
import numpy as np


_NAME_BY_INDEX = ("small", "medium", "large")


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
                 vlim: float | None) -> None:
    """Render the y-midplane (x-z slice) as a PNG with a diverging colormap.

    The strain field is signed ([-alpha, alpha]); a diverging colormap
    centered at zero reads better than a sequential one. The y-midplane
    is chosen so the dilation axis (when --axis x) varies horizontally
    and the orthogonal direction varies vertically.
    """
    Ny = field_3d.shape[1]
    slice_2d = field_3d[:, Ny // 2, :]
    if vlim is None:
        vlim = max(1e-6, float(np.abs(slice_2d).max()))
    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    im = ax.imshow(slice_2d, origin="lower", aspect="equal",
                   cmap="RdBu_r", vmin=-vlim, vmax=vlim)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$\varepsilon_{\mathrm{axial}}$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-dir", type=Path, required=True,
                        help="Directory containing the input .vtk files.")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Directory to write the PNGs into.")
    parser.add_argument("--files", nargs="+", required=True,
                        help="Input .vtk basenames in order (small, medium, large).")
    args = parser.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # First pass: determine a shared color limit so the three PNGs are
    # visually comparable.
    fields = []
    for name in args.files:
        _, data = parse_vtk(args.in_dir / name)
        fields.append(data)
    vlim = max(float(np.abs(f).max()) for f in fields)

    for idx, (name, data) in enumerate(zip(args.files, fields)):
        out_name = (_NAME_BY_INDEX[idx] if idx < len(_NAME_BY_INDEX)
                    else f"snapshot_{Path(name).stem}") + ".png"
        out_path = args.out_dir / out_name
        title = f"Axial strain ({name})"
        render_slice(data, title, out_path, vlim)
        print(f"wrote {out_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
