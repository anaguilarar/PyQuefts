# Step 3 — PyQUEFTS architecture

## Design principles

1. **Core is pure NumPy and array-shaped.** Every scientific function accepts scalars or broadcastable arrays and returns arrays of the broadcast shape. The same code therefore serves a single point, a table of points and a whole raster (after flattening or via `xarray.apply_ufunc`). No Python loops over pixels.
2. **Scientific core knows nothing about files, units or geography.** Unit conversion, depth aggregation and georeferencing live in `pyquefts.geo`; the core receives values in model units (SOC g/kg, P mg/kg, Kex mmol/kg, pH, °C, kg/ha).
3. **Parameters are explicit, immutable dataclasses** with the exact R-QUEFTS names, so R and Python configurations can be compared field by field.
4. **Intermediate results are first-class outputs.** `run()` returns a `QueftsResult` holding supplies, pairwise uptakes, uptake, yield ranges, six pairwise yields, limited yields and gaps.
5. **Extensions are separated from the reference model.** Anything that is not in R-QUEFTS lives in `pyquefts.analysis` or `pyquefts.geo` and is documented as an extension.

## Package layout

```
pyquefts/
  __init__.py            public API re-exports
  core/
    parameters.py        SoilParameters, CropParameters, Biomass, Fertilizer, UPTAKE_ADJUST_DEFAULT
    supply.py            nut_supply1 (Janssen 1990), nut_supply2 (Sattari 2014), pH factor tables
    conversions.py       Mehlich-3 -> Olsen P (McDowell et al. 2023); registry of P methods
    model.py             uptake(), pair_yield(), requirements(), run_quefts()  — vectorized core
    result.py            QueftsResult container (dict-like, to_dict / to_dataframe / to_dataset)
  crops.py               crop table access (get_crop, list_crops), leaf/stem ratio helpers
  point.py               Module 1: Quefts class (R-like object) and run_batch()
  fertilizer.py          fertilizers table, nutrient_rates, fert_app (LP), opt_app (profit optimum)
  reverse.py             rev_supply (fit base supplies to omission-trial data)
  analysis.py            Extensions: response curves/surfaces, uptake–yield envelopes, balanced requirements
  plot.py                Extensions: matplotlib plots of response curves/surfaces, envelopes (optional dependency)
  geo/
    cube.py              SoilCube: open NetCDF/Zarr/GeoTIFF stack, variable mapping, unit conversion,
                         depth aggregation to a target layer (default 0–20 cm)
    extract.py           Module 2: extract_points(cube, lon, lat) -> DataFrame; run_points(...); response_points(...)
    spatial.py           Module 3: run_spatial(cube, crop, ...) -> xarray.Dataset (dask-chunked);
                         response_spatial (rate dimension); write_netcdf, write_geotiff
    climate.py           growing-season mean temperature from the CHIRTS-ERA5 daily Zarr (Hugging Face),
                         regridding onto the soil grid, point sampling
  data/
    quefts_crop_pars.csv     verbatim from Rquefts
    fertilizers.csv          verbatim from Rquefts
tests/                   pytest unit tests (core, supply, crops, geo on synthetic cubes)
validation/              R script producing reference outputs + Python comparison script + report
examples/                runnable scripts for the three modes
docs/                    this documentation
```

## Data flow of the three modes

```
Mode 1 (point):   user values ──► nut_supply1/2 ──► SoilParameters ──► run_quefts ──► QueftsResult
Mode 2 (coords):  cube ──► nearest cell ──► depth-aggregate ──► unit convert ──► (same as mode 1, vectorized over points)
Mode 3 (spatial): cube[bbox/mask] ──► depth-aggregate (xarray) ──► unit convert ──► apply_ufunc(run_quefts, dask) ──► Dataset ──► NetCDF/GeoTIFF
```

## Key interfaces

```python
from pyquefts import Quefts, get_crop, nut_supply1, nut_supply2

# Mode 1
sup = nut_supply2(temp=25, ph=6.0, soc=15, kex=1.0, polsen=10, ptotal=500)
q = Quefts(soil=sup.to_soil(), crop=get_crop("maize"),
           fert=dict(N=60, P=10, K=60), biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
res = q.run()           # QueftsResult; res["store_lim"], res.to_dict()

# Mode 2
from pyquefts.geo import SoilCube, run_points
cube = SoilCube.open("hnd_sh_soilgrids_format.nc")          # auto-detects variables/units via a mapping
df = run_points(cube, lon=[-87.6], lat=[14.45], crop="maize", temp=24, store_att=6000, ...)

# Mode 3
from pyquefts.geo import run_spatial
ds = run_spatial(cube, crop="maize", temp=24, store_att=6000, bbox=(-88, 14, -87, 15), fert=dict(N=0,P=0,K=0))
ds.to_netcdf("quefts_maize.nc")
```

## Vectorization notes

* `uptake()` evaluates the three situations (A, B, C) with `numpy.where`; divisions that are undefined in a masked-out branch are computed under `numpy.errstate(all="ignore")`.
* `pair_yield()` implements the three-mask reduction of the C++ control flow (see `02_r_reference_analysis.md`, section 4), including the re-computation of the upper limit when the modelled nutrient is itself the most limiting.
* Organ selection (`store_att > 200` etc.) is also a mask so that attainable-yield rasters that contain tiny values do not break the vectorization.
* Spatial runs use `xarray.apply_ufunc(..., dask="parallelized")` on 2-D chunks; each chunk flattens to 1-D, runs the core once and reshapes. Memory stays bounded by chunk size (the full Honduras cube is 64 GB on disk).

## Units and the Honduras cube

The provided cube (`hnd_sh_soilgrids_format.nc`) stores SoilGrids-style integer-scaled units. `SoilCube` carries a mapping table (`VARIABLE_MAP`) that converts to model units and that can be overridden per dataset:

| cube variable | cube unit | model variable | conversion |
|---|---|---|---|
| `phh2o` | pH × 10 | pH | ÷ 10 |
| `soc` | dg/kg | SOC g/kg | ÷ 10 |
| `nitrogen` | cg/kg | Norg g/kg | ÷ 100 |
| `k` | cmol(+)/kg | Kex mmol/kg | × 10 |
| `p` | mg/kg Mehlich-3 | P-Olsen mg/kg | Olsen = 0.47 × M3 + 2.4 (McDowell et al. 2023, non-calcareous) |

Depth layers `0-5`, `5-15`, `15-30`, `30-60` cm are aggregated to the 0–20 cm reference layer with thickness weights (5, 10, 5) by default. The `p` and `k` layers are the same at every depth (topsoil value repeated, per the file metadata), so aggregation is a no-op for them.

**Note on P.** The cube's `p` variable is Mehlich-3 extractable P from the Suelos de Honduras survey (confirmed by the data owner; the file metadata does not say so). QUEFTS Step 1 is calibrated on P-Olsen, so `SoilCube.open(..., p_method="mehlich3")` converts it with the McDowell et al. (2023) non-calcareous regression (`pyquefts.core.conversions`); the default `p_method="olsen"` leaves values untouched. A QUEFTS-specific transfer function (Breure et al. 2022) also needs Mehlich-3 Al and cannot be applied to this cube. Total P is not in the cube, so `nut_supply2` needs a `ptotal_ratio` (default 55 × P-Olsen following the Rquefts documentation; Sattari's 95 is available); the ratio is applied to the converted Olsen value.

## Temperature input

The cube has no temperature. `pyquefts.geo.climate.season_mean_temperature` computes the growing-season mean from the CHIRTS-ERA5 daily Tmax/Tmin Zarr on Hugging Face (0.05°), as a multi-year climatology or for a given season, and `run_spatial`/`run_points` accept the resulting `(lat, lon)` field directly (bilinear regridding / nearest sampling). See `docs/06_temperature_from_chirts.md`.
