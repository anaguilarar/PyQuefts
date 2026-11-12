"""Step 6: performance benchmarks of the PyQUEFTS core and the spatial module.

Usage: python benchmark.py [--cube PATH] [--out results.md]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyquefts as pq  # noqa: E402

DEFAULT_CUBE = r"D:/OneDrive - CGIAR/projects/spatial_files/raster/soil/hnd_sh_soilgrids_format.nc"


def bench_core(sizes=(1, 1_000, 100_000, 1_000_000, 4_000_000)) -> list[dict]:
    rng = np.random.default_rng(0)
    rows = []
    crop = pq.get_crop("maize")
    for n in sizes:
        soil = pq.SoilParameters(N_base_supply=rng.uniform(0, 200, n), P_base_supply=rng.uniform(0, 40, n),
                                 K_base_supply=rng.uniform(0, 300, n))
        fert = pq.Fertilizer(N=rng.uniform(0, 150, n), P=rng.uniform(0, 40, n), K=rng.uniform(0, 100, n))
        biom = pq.Biomass.from_yield(rng.uniform(1000, 12000, n), 0.46, 0.56)
        for keep in (False, True):
            t = time.perf_counter()
            reps = 5 if n <= 1000 else 1
            for _ in range(reps):
                pq.run_quefts(soil, crop, fert, biom, keep_intermediate=keep)
            dt = (time.perf_counter() - t) / reps
            rows.append({"n_elements": n, "keep_intermediate": keep, "seconds": dt,
                         "elements_per_s": n / dt if dt > 0 else np.inf})
    return rows


def bench_python_loop(n=2000) -> float:
    """Per-call overhead of the object interface (equivalent of an R for-loop over run())."""
    q = pq.Quefts()
    t = time.perf_counter()
    for i in range(n):
        q.set(N=float(i % 200))
        q.run(keep_intermediate=False)
    return (time.perf_counter() - t) / n


def bench_spatial(cube_path: str, bboxes, chunk_sizes=(512, 1024, 2048)) -> list[dict]:
    from pyquefts.geo import SoilCube, run_spatial

    rows = []
    for cs in chunk_sizes:
        cube = SoilCube.open(cube_path, chunks={"depth": -1, "y": cs, "x": cs})
        for bbox in bboxes:
            if cs != 1024 and bbox is bboxes[-1]:
                continue  # the largest region only with the default chunk size
            t = time.perf_counter()
            ds = run_spatial(cube, "maize", store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                             fert=dict(N=0, P=0, K=0), temp=24.0, bbox=bbox, keep_inputs=False, compute=True)
            dt = time.perf_counter() - t
            ncell = int(ds.store_lim.size)
            rows.append({"chunk": cs, "bbox": bbox, "cells": ncell, "seconds": dt, "cells_per_s": ncell / dt,
                         "valid_frac": float(ds.store_lim.notnull().mean())})
            print(rows[-1], flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cube", default=DEFAULT_CUBE)
    ap.add_argument("--out", default=None)
    ap.add_argument("--skip-spatial", action="store_true")
    a = ap.parse_args()

    core = bench_core()
    loop = bench_python_loop()
    lines = ["# PyQUEFTS performance benchmark", "", f"NumPy {np.__version__}", "", "## Core (vectorized run_quefts)", "",
             "| n elements | intermediates | seconds | elements/s |", "|---|---|---|---|"]
    for r in core:
        lines.append(f"| {r['n_elements']:,} | {r['keep_intermediate']} | {r['seconds']:.4f} | {r['elements_per_s']:,.0f} |")
    lines += ["", f"Object interface, one call per scenario (Python loop): {loop*1e3:.3f} ms per run", ""]
    print("\n".join(lines))

    if not a.skip_spatial and Path(a.cube).exists():
        bboxes = [(-87.70, 14.40, -87.60, 14.50), (-88.0, 14.0, -87.5, 14.5), (-88.5, 13.5, -87.5, 14.5)]
        sp = bench_spatial(a.cube, bboxes)
        lines += ["## Spatial (Honduras cube, Sattari supply, maize, dask chunks)", "",
                  "| chunk | bbox | cells | seconds | cells/s | valid frac |", "|---|---|---|---|---|---|"]
        for r in sp:
            lines.append(f"| {r['chunk']} | {r['bbox']} | {r['cells']:,} | {r['seconds']:.1f} | {r['cells_per_s']:,.0f} | {r['valid_frac']:.2f} |")
    if a.out:
        Path(a.out).write_text("\n".join(lines), encoding="utf-8")
        print("written", a.out)


if __name__ == "__main__":
    main()
