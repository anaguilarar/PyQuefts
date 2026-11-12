# Growing-season temperature from CHIRTS-ERA5

## Why temperature is needed

The Sattari et al. (2014) recalibration used by `nut_supply2` replaces the constant N-mineralization coefficient 6.8 of Janssen (1990) by a temperature-dependent one, `alphaN = 2 * 2^((T - 9)/9)` (see `docs/01_scientific_review.md`, section 2.2), where T is the **mean air temperature of the growing season** (°C). The Honduras soil cube has no temperature layer, so T must come from a climate dataset. At 24 °C the coefficient is 6.35; at 28 °C it is 8.6, so a 4 °C difference changes the native N supply by about 35 %: the choice matters.

## Source

CHIRTS-ERA5 daily maximum and minimum 2 m air temperature, Climate Hazards Center (UC Santa Barbara), redistributed as a Zarr v3 datacube on Hugging Face by A. Aguilar (CC-BY-4.0):

* https://huggingface.co/datasets/aaguilar90/chirts-era5-daily
* 0.05° global grid (60 S-70 N), daily 1981-01-01 to present, variables `Tmax` and `Tmin` in °C (stored as int16 in 0.1 °C; xarray applies the scale factor).
* Chunks (128 days, 300 lat, 300 lon); only the chunks intersecting the query are transferred. A whole-Honduras season loads in about 2 s.
* Caveat from the dataset card: this is CHIRTS-ERA5 (experimental tree), not the station-blended CHIRTS-daily v1.0 (1983-2016); do not mix the two.

## What PyQUEFTS computes

`pyquefts.geo.climate.season_mean_temperature(start, end, years, bbox)`:

1. selects the season window `start..end` (``"MM-DD"``) in each requested year (a window ending before it starts is taken to cross New Year, e.g. `11-01..01-31`);
2. computes the daily mean as `(Tmax + Tmin) / 2`;
3. averages over all days of all selected years, giving a climatological growing-season mean (default years 1991-2020; pass a single year for a specific season);
4. fills water / nodata cells with the nearest valid cell so that coastal soil pixels get a value (`fill=True`);
5. returns a `(lat, lon)` DataArray in °C with the window and years in its attributes.

`run_spatial(..., temp=<that DataArray>)` interpolates it bilinearly onto the soil grid (`regrid_to`, 0.05° to 0.00028°) and keeps it chunked like the cube; `run_points(..., temp=<DataArray>)` samples the nearest climate cell per point. The temperature used is stored in the output (`temp` variable / column).

Using `(Tmax + Tmin)/2` as the daily mean is the standard approximation when hourly data are unavailable; Sattari et al. (2014) used mean growing-season temperature from station or reanalysis data without specifying the daily aggregation.

## Honduras defaults used in the examples

| Season | Window | Why |
|---|---|---|
| Primera (maize, beans) | 05-15 to 08-31 | main rainy-season sowing in May, harvest August |
| Postrera | 09-01 to 12-15 | second season |

Users should set the window to the crop calendar of the region and crop; the `years` argument allows a specific season (for example `years=[2024]`) or a climatology.

For the Comayagua test window (14.4-14.5 N, 87.6-87.7 W) the primera 1991-2020 mean is about 25-27 °C in the valley; CHIRTS resolves the cooler mountains around it (about 20 °C), which the previous constant of 24 °C ignored.

## How to use

```python
from pyquefts.geo import SoilCube, season_mean_temperature, run_spatial, run_points

temp = season_mean_temperature(start="05-15", end="08-31", years=range(1991, 2021),
                               bbox=(-89.4, 12.9, -83.1, 16.6))      # whole Honduras, ~2 s
cube = SoilCube.open("hnd_sh_soilgrids_format.nc", chunks={"depth": -1, "y": 1024, "x": 1024},
                     p_method="mehlich3")
ds = run_spatial(cube, "maize", store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                 fert=dict(N=0, P=0, K=0), temp=temp, bbox=(-88, 14, -87.5, 14.5), compute=True)
df = run_points(cube, lon=[-87.65], lat=[14.45], crop="maize", temp=temp, store_att=6000,
                leaf_ratio=0.46, stem_ratio=0.56)
```
A local copy of the store, or any other `(time, lat, lon)` dataset with daily `Tmax`/`Tmin` in °C, can be passed as the first argument of `season_mean_temperature` (variable and dimension names are parameters).

## Tests

`tests/test_climate.py` checks the season selection (including windows crossing New Year), the climatology averaging, nodata filling, bilinear regridding onto a dask-chunked soil grid, point sampling, and the end-to-end use in `run_spatial`/`run_points` on synthetic data. One test queries the Hugging Face store for a Honduras subset and runs only when `PYQUEFTS_NETWORK_TESTS=1` is set.
