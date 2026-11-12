# Step 4-6 — Validation against R-QUEFTS and performance

## 1. Validation framework

`validation/` contains a reproducible two-sided framework:

| File | Purpose |
|---|---|
| `generate_r_reference.R` | Runs the installed `Rquefts` package (CRAN 1.2-8, R 4.4.1) on randomized and documented scenarios and writes CSVs: `r_reference_model.csv` (5000 model runs), `r_reference_supply.csv` (2000 soils for `nutSupply1/2`), `r_reference_examples.csv` (tutorial / manual examples). |
| `compare_with_r.py` | Runs PyQUEFTS on the same inputs (vectorized, one call per crop) and compares every output at `rtol = 1e-9, atol = 1e-7`. Writes `VALIDATION_REPORT.md`. |
| `benchmark.py` | Step 6 timing of the core and of spatial runs on the Honduras cube. Writes `BENCHMARK.md`. |
| `tests/` | 41 pytest unit tests: documented R outputs, the Janssen (1990) Table VI worked example, vectorized-vs-scalar consistency, organ branches, NaN handling, supply functions, geospatial extraction and spatial runs on a synthetic cube, fertilizer/reverse/analysis tools. |

Scenario design for the randomized runs (seed 20251112):

* all 28 crop parameter sets plus the R default set;
* base supplies N 0-200, P 0-40, K 0-300 kg/ha with 10 % zeros and 5 % values up to 5x the range;
* recoveries N 0.2-0.8, P 0.05-0.4, K 0.2-0.8;
* fertilizer N 0-250, P 0-80, K 0-200 kg/ha with 30 % zeros;
* attainable storage yield 500-15000 kg/ha; 8 % of runs with storage organ < 200 kg/ha (half of those also with stems < 200) to exercise the leaf- and stem-modelled branches; leaf and stem biomass 0.1-1.2x the storage yield;
* season length in {60, 90, 110, 120, 150, 180, 240, 400} days;
* supply soils: pH 3.5-9.5 (outside the model's validity range on purpose, to test clamping), SOC 0-80 g/kg, Kex 0-40 mmol/kg, P-Olsen 0-40 mg/kg, total P 0-2000 mg/kg, T 10-35 C.

## 2. Results

From `validation/VALIDATION_REPORT.md` (generated 2025-11-12):

| Comparison | n | mismatches | max relative diff |
|---|---|---|---|
| 12 model outputs (supply, uptake, leaf/stem/store yield, gaps) | 5000 | 0 | 5e-15 (2e-13 for P, values near zero) |
| `nut_supply1` with and without total P | 2000 | 0 | 5e-12 |
| `nut_supply2` | 2000 | 0 | 1e-13 |
| Documented examples (default run, barley, rice N series) | 7 cases x 12 outputs | 0 | abs 5e-12 |

The differences are floating-point rounding (different operation order in NumPy vs C++). The intermediate quantities of each step are also validated indirectly: the Janssen (1990) Table VI example is reproduced step by step in `tests/test_core.py::test_janssen_1990_table_vi_soil_a` (supplies 102/9.5/125 → pairwise uptakes 94, 101, 9.3, 113, 108 → yield ranges → six pairwise yields → YE 4279 kg/ha, all within rounding of the published table).

### Optimization-based functions

| Function | R result | PyQUEFTS result | Explanation |
|---|---|---|---|
| `fert_app` (LP, manual example) | urea 86.43269, TSP 123.51847, K-mag 0, NPK 301.20482 (and the four other targets) | identical to 5 decimals | both solve the same LP (limSolve::linp vs scipy HiGHS) |
| `opt_app`, dm value 0.75 | N 211.38, P 13.36, K 25.40; urea 393, NPK 153; NR 304 | identical | same optimum found |
| `opt_app`, dm value 0.25 | urea 28 kg, NR 2 | no fertilizer, NR 0 | flat objective near zero; R `optim(L-BFGS-B, ndeps=1, fnscale=-0.1)` and SciPy L-BFGS-B stop at different points of a plateau. Both are economically equivalent (NR within 2 units). |
| `rev_supply` (potato trial example) | 35.23, 8.01, 385.07 | 35.23, 8.00, 385.63 | same minimum of a flat SSE surface in K (the trial has no K response); differences are optimizer tolerance. A synthetic recovery test with known supplies passes at 5 %. |

### Documented behavioural differences

1. `season_factor` interpolates the full `UptakeAdjust` table; the C++ `approx()` drops the last row. Identical for the default table (verified above); for a custom table whose last two factors differ, R would ignore the last breakpoint.
2. `nut_supply2(ptotal=None)` estimates total P as `ratio * P-Olsen` (default 55, Rquefts documentation; 95 per Sattari et al. 2014). R requires total P.
3. PyQUEFTS propagates NaN from any supply or biomass input to all outputs; R `run()` would return numbers computed from NaN comparisons, and R `predict()` only skips cells with NaN N supply.
4. The `leaf_lim` fallback for zero yield reproduces the R expression `(NminVeg + 2*NminVeg)` verbatim (see `docs/01_scientific_review.md`, section 9).

## 3. Edge cases and input validation covered by tests

* zero supplies and zero fertilizer (yield 0, finite partition);
* NaN inputs (per element);
* storage organ or stems below 200 kg/ha (organ branches);
* pH outside 4.5-7 (clamping in `nut_supply1`, table clamping in `nut_supply2`);
* decreasing `UptakeAdjust` table (error), unknown crop (error), unknown units in a cube (error);
* cube depth labels not covering the target layer (error); coordinates outside the cube (NaN row);
* custom unit transforms and the Janssen-vs-Sattari clamping behaviour on grids.

## 4. Geospatial consistency

`tests/test_geo.py` builds a synthetic SoilGrids-style cube (integer-scaled units, three depths, NaN cells, descending latitude) and checks that:

* the thickness-weighted 0-20 cm aggregation equals the hand-computed weights (0.25, 0.5, 0.25 for 0-5, 5-15, 15-30 cm);
* `extract_points` returns exactly the aggregated value of the nearest cell and NaN outside the grid;
* `run_points` and `run_spatial` give the same numbers as the point mode for the same cell;
* a raster-valued attainable yield and a bounding box are honoured;
* NetCDF and 12-band GeoTIFF round-trips preserve values and CRS.

On the real Honduras cube (`hnd_sh_soilgrids_format.nc`), two test points and a 0.1 x 0.1 degree subset (131k cells) were run; the subset run takes about 4 s after the lazy open, and 12 % of cells are NaN (no-data in the source).

## 5. Performance (Step 6)

Measured on this machine (Windows 11, NumPy 2.4, single process; full table in `validation/BENCHMARK.md`):

| Workload | Result |
|---|---|
| Core, 1 element | 6 ms per call (fixed overhead of broadcasting ~40 arrays) |
| Core, 100 000 elements | 0.2 s (about 500 000 elements/s) |
| Core, 1 000 000 elements | 3.3 s without intermediates, 4.7 s with |
| Core, 4 000 000 elements | 14.5 s (about 275 000 elements/s, memory-bandwidth bound) |
| Object interface in a Python loop | 1.7 ms per `Quefts.run()` (R-style loop) |
| Point extraction, 3 points (Module 2) | 0.6 s without dask, 2.4 s with dask chunks (one storage chunk read per point and variable) |
| Point extraction, 500 points scattered over Honduras | 134 s (about 1000 compressed-chunk reads at ~0.13 s each; dense point sets are much cheaper) |
| Spatial, 0.1 x 0.1 deg (131 k cells), chunk 1024 | 0.75 s |
| Spatial, 0.5 x 0.5 deg (3.3 M cells), chunk 512 / 1024 / 2048 | 7.7 / 8.1 / 10.7 s (about 400 000 cells/s) |
| Spatial, 1 x 1 deg (13.1 M cells), chunk 1024 | 31 s (422 000 cells/s) |

Extrapolated, the full Honduras cube (287 M cells, of which roughly half are land) would take about 10-15 minutes in a single process; dask's threaded scheduler parallelizes chunks automatically when several cores are available. Guidance:

* The core is fully vectorized; cost is linear in the number of elements with no Python loop. Keeping the ~40 intermediate arrays roughly doubles memory; use `keep_intermediate=False` (the default in the spatial module when only R outputs are requested) for large grids.
* The object interface (`Quefts.run()` in a Python loop) has a fixed overhead per call; for sweeps, pass arrays instead of looping.
* Spatial runs are chunked with dask: each (y, x) chunk is flattened and run once. Chunk size trades scheduler overhead (small chunks) against memory (large chunks); 1024 x 1024 is a good default for the Honduras cube (chunks on disk are 1024 x 1024 x 1 depth). Reading and decompressing the NetCDF is the dominant cost, not the model.
* For point queries open the cube without dask (`chunks=None`); the extractor groups points by storage chunk and reads each chunk once. For thousands of scattered points a spatial run over the bounding region followed by sampling can be faster.
* For whole-country runs write the result lazily (`compute=False` then `write_netcdf`) so that dask streams chunks to disk; or convert the source cube to Zarr once for faster chunked reads.
