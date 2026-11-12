# PyQUEFTS

Python implementation of **QUEFTS** (Quantitative Evaluation of the Fertility of Tropical Soils; Janssen et al., 1990), ported from and validated against the R package [`Rquefts`](https://github.com/cropmodels/Rquefts) (v1.2-8, Hijmans et al.), with vectorized NumPy calculations and geospatial support (NetCDF / Zarr soil data cubes, GeoTIFF output).

QUEFTS (1) estimates the native N, P and K supply of a soil from a few topsoil properties and (2) computes nutrient uptake, nutrient-limited yield and fertilizer requirements from that supply, crop parameters, fertilizer applied and the attainable yield.

## Status

| Component | State |
|---|---|
| Core model (steps 2-4), crop table, season adjustment | Complete. Bit-for-bit agreement with Rquefts on 5000 randomized scenarios (max relative diff 5e-15). See `validation/VALIDATION_REPORT.md`. |
| Supply functions `nut_supply1` (Janssen 1990) and `nut_supply2` (Sattari 2014) | Complete, validated on 2000 random soils. |
| Fertilizer tools (`fert_app`, `opt_app`, `nutrient_rates`), reverse QUEFTS (`rev_supply`) | Complete. `fert_app` identical to R; the two optimizers agree to optimizer tolerance (see validation doc). |
| Module 1: point simulation | Complete (`Quefts`, `run_point`, `run_batch`). |
| Module 2: coordinates from a soil cube | Complete (`pyquefts.geo.run_points`). |
| Module 3: spatial simulation | Complete (`pyquefts.geo.run_spatial`, dask-chunked, NetCDF/GeoTIFF). |
| Extensions (response curves/surfaces, uptake-yield envelopes, balanced uptake) | Complete, in `pyquefts.analysis`; documented as extensions. |
| Plotting (`result.plot()`, `Quefts.plot()`, per-point and regional responses from a cube, surfaces, envelopes) | Complete, in `pyquefts.plot`, `geo.response_points`, `geo.response_spatial` (needs `pip install pyquefts[plot]`). See `docs/07_plotting.md`. |
| Input preparation: Mehlich-3 to Olsen P (McDowell et al. 2023), growing-season temperature from CHIRTS-ERA5 | Complete, in `pyquefts.core.conversions` and `pyquefts.geo.climate`. |

## Installation

Requires Python >= 3.9. The core (point simulations) depends only on NumPy, pandas and SciPy. The remote CHIRTS-ERA5 temperature store is Zarr v3 and needs `zarr >= 3` (Python >= 3.11); everything else, including local NetCDF/GeoTIFF cubes, works on 3.9.

Install from GitHub (https://github.com/anaguilarar/PyQuefts):

```bash
pip install "git+https://github.com/anaguilarar/PyQuefts.git"            # core only
pip install "pyquefts[all] @ git+https://github.com/anaguilarar/PyQuefts.git"   # with geo, plot and test extras
```

Or clone it for development (editable install):

```bash
git clone https://github.com/anaguilarar/PyQuefts.git
cd PyQuefts
pip install -e ".[all]"
```

Optional extras:

| Extra | Adds | Needed for |
|---|---|---|
| `geo` | xarray, dask, netCDF4, rioxarray, rasterio, zarr | soil cubes (`pyquefts.geo`), spatial runs, CHIRTS temperature |
| `plot` | matplotlib | `result.plot()`, `pyquefts.plot` |
| `test` | pytest, tabulate | running the test suite |
| `all` | all of the above | |

Check the installation:

```bash
python -c "import pyquefts as pq; print(pq.__version__, pq.Quefts(crop='barley').run())"
pytest            # from the cloned folder; set PYQUEFTS_NETWORK_TESTS=1 to also query the CHIRTS store
```

The R reference implementation is only needed to regenerate the validation data (`validation/generate_r_reference.R`, Rquefts 1.2-8).

## Quick start

```python
import pyquefts as pq

# Step 1: native supply (kg/ha, 120-day season) from 0-20 cm topsoil properties
sup = pq.nut_supply2(temp=25, ph=6.0, soc=15, kex=1.0, polsen=10, ptotal=500)

# Steps 2-4: nutrient-limited yield for barley with 100 kg N/ha
q = pq.Quefts(soil=sup.to_soil(), crop="barley", fert=dict(N=100, P=0, K=0),
              biom=dict(leaf_att=2200, stem_att=2700, store_att=4800, SeasonLength=110))
res = q.run()
res["store_lim"], res["N_uptake"], res["P_gap"]      # grain yield, N uptake, P fertilizer needed
res.to_dict(r_outputs_only=True)                     # the 12 outputs of R `run()`
res["YND"], res["Y_NP"]                              # intermediates: yield ranges, pairwise yields
```

All inputs may be arrays, so a fertilizer sweep is a single call:

```python
import numpy as np
res = pq.Quefts(crop="rice", fert=dict(N=np.arange(0, 201, 10), P=0, K=0),
                biom=dict(leaf_att=2651, stem_att=5053, store_att=8208)).run()
res["store_lim"]   # one yield per N rate
```

### Coordinates from a soil cube (Module 2)

```python
from pyquefts.geo import SoilCube, run_points
cube = SoilCube.open("hnd_sh_soilgrids_format.nc", chunks={"depth": -1, "y": 1024, "x": 1024},
                     p_method="mehlich3")   # the Honduras survey P is Mehlich-3
print(cube)   # shows how variables/units were mapped
df = run_points(cube, lon=[-87.65, -86.2], lat=[14.45, 14.8], crop="maize", temp=24,
                store_att=6000, leaf_ratio=0.46, stem_ratio=0.56, fert=dict(N=0, P=0, K=0))
```

### Spatial simulation (Module 3)

```python
from pyquefts.geo import run_spatial, season_mean_temperature, write_netcdf, write_geotiff
temp = season_mean_temperature(start="05-15", end="08-31", years=range(1991, 2021),
                               bbox=(-88, 14, -87.5, 14.5))     # CHIRTS-ERA5 climatology, 0.05 deg
ds = run_spatial(cube, "maize", store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                 fert=dict(N=60, P=15, K=30), temp=temp, bbox=(-88, 14, -87.5, 14.5), compute=True)
write_netcdf(ds, "quefts_maize.nc")
write_geotiff(ds, "quefts_maize.tif", variables=["store_lim", "N_gap", "P_gap", "K_gap"])
```
`store_att`, `temp`, fertilizer rates and season length may each be scalars or rasters; a field on another grid (such as the 0.05° CHIRTS temperature) is interpolated onto the soil grid automatically.

## Units

| Input | Unit | Note |
|---|---|---|
| pH | pH (H2O) | |
| SOC | g/kg | organic carbon, 0-20 cm |
| Kex | mmol/kg | exchangeable K (1 cmol(+)/kg = 10 mmol/kg) |
| P-Olsen, total P | mg/kg | Mehlich-3 P can be converted with `p_method="mehlich3"` (McDowell et al. 2023) |
| temperature | C | growing-season mean (Sattari 2014 supply only); derive it from CHIRTS-ERA5 with `season_mean_temperature` |
| supplies, uptake, fertilizer, gaps | kg/ha | fertilizer as elemental N, P, K |
| biomass / yield | kg/ha dry matter | leaves, stems, storage organ |

`SoilCube` converts SoilGrids-style storage units (pH*10, dg/kg, cg/kg, cmol(+)/kg) automatically from the `units` attributes and lets you override any mapping, including a callable transform (for example a Mehlich-3 to Olsen conversion).

## Documentation

* `docs/01_scientific_review.md` — equations, parameters and calculation workflow from Janssen (1990), Smaling & Janssen (1993), Sattari et al. (2014).
* `docs/02_r_reference_analysis.md` — structure of the R implementation and the mapping to PyQUEFTS.
* `docs/03_architecture.md` — package design, the three input modes, vectorization notes, units of the Honduras cube.
* `docs/04_validation_and_performance.md` — validation framework, results, explained differences, benchmarks.
* `docs/05_mehlich3_to_olsen_conversion.md` — why and how Mehlich-3 P is converted to Olsen P, with sources.
* `docs/06_temperature_from_chirts.md` — growing-season temperature from the CHIRTS-ERA5 daily cube.
* `docs/07_plotting.md` — response-curve, response-surface and uptake-yield plots.
* `docs/08_prompt_dssat_quefts_coupling.md` — prompt for coupling DSSAT (WeatherSoilDataProcessor) yield ensembles to PyQUEFTS for probabilistic recommendations.
* `examples/` — runnable scripts for the three modes.

## Differences from R-QUEFTS (all documented)

* Season-length interpolation uses the full `UptakeAdjust` table (the C++ `approx()` drops the last row; identical for the default table).
* `nut_supply2` can estimate total P from P-Olsen with a configurable ratio when total P is absent (R requires total P).
* `pyquefts.analysis`, the depth-aggregation / unit-conversion layer and the Mehlich-3 to Olsen P conversion (`pyquefts.core.conversions`, McDowell et al. 2023) are additions; they do not change the model equations.

## License

GPL-3.0-or-later, as the reference implementation (see `docs/RQUEFTS_LICENSE`). Crop and fertilizer tables are copied verbatim from Rquefts.

## References

* Janssen, B.H., Guiking, F.C.T., van der Eijk, D., Smaling, E.M.A., Wolf, J., van Reuler, H., 1990. A system for quantitative evaluation of the fertility of tropical soils (QUEFTS). Geoderma 46, 299-318.
* Smaling, E.M.A., Janssen, B.H., 1993. Calibration of QUEFTS, a model predicting nutrient uptake and yields from chemical soil fertility indices. Geoderma 59, 21-44.
* Sattari, S.Z., van Ittersum, M.K., Bouwman, A.F., Smit, A.L., Janssen, B.H., 2014. Crop yield response to soil fertility and N, P, K inputs in different environments: Testing and improving the QUEFTS model. Field Crops Research 157, 35-46.
* McDowell, R.W., Pletnyakov, P., Noble, A., Haygarth, P.M., 2023. A global database of soil plant available phosphorus. Scientific Data 10, 125. (Mehlich-3 to Olsen P conversion)
* Breure, T.S., Van Eynde, E., Kempen, B., Comans, R.N.J., Hoffland, E., 2022. Transfer functions for phosphorus and potassium soil tests and implications for the QUEFTS model. Geoderma 406, 115458.
* Funk, C. et al. CHIRTS-daily / CHIRTS-ERA5, Climate Hazards Center, UC Santa Barbara; Zarr redistribution: https://huggingface.co/datasets/aaguilar90/chirts-era5-daily (CC-BY-4.0).
* Hijmans, R.J., Pypers, P., Wolff, J. Rquefts: Quantitative Evaluation of the Native Fertility of Tropical Soils. R package version 1.2-8. https://CRAN.R-project.org/package=Rquefts
