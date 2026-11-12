# Plotting nutrient responses

`pyquefts.plot` is an extension: it only draws outputs of the unchanged model. It needs matplotlib (`pip install pyquefts[plot]`), which is imported lazily, so the rest of the package works without it.

## The Rquefts loop, vectorized

The Rquefts tutorial builds a response curve with a loop:

```r
for (i in 1:length(N)) { q["N"] <- N[i]; results[i,] <- run(q) }
plot(N, results[, "store_lim"])
```

In PyQUEFTS every parameter may be an array, so one run gives the whole curve, and the inputs are recorded on the result so that `plot()` knows what varied:

```python
import numpy as np
from pyquefts import Quefts, nut_supply2

sup = nut_supply2(temp=26, ph=6.0, soc=15, kex=1.0, polsen=10)
q = Quefts(soil=sup.to_soil(), crop="maize",
           fert=dict(N=np.arange(0, 160, 10), P=10, K=60),
           biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
res = q.run()                      # 16 elements
res["store_lim"]                   # yield per N rate
res.to_dataframe(include_inputs=True)   # N column + outputs
res.plot()                         # store_lim vs N
res.plot(y=["store_lim", "N_uptake", "N_gap"])
```

`res.plot()` finds the single input that varies (`N`, `P`, `K`, `store_att`, `N_base_supply`, a recovery, ...). When several vary, name them: `res.plot(x="N", hue="store_att")` draws one line per attainable yield. `x` may also be an output, e.g. `res.plot(x="N_uptake", y="store_lim")`.

Fertilizer keys are case-insensitive (`dict(N=..., P=10, k=60)` works).

## Three panels at once

```python
q.plot(rates=np.arange(0, 201, 10))          # N, P and K panels; other two nutrients held at q.fert
q.response("P", [0, 10, 20, 40])             # the same as a DataFrame (analysis.response_curve)
```

`plot_response_curves(crop=..., soil=..., fert=..., biom=..., rates={"N": ..., "P": ..., "K": ...})` does the same without a `Quefts` object.

## Surfaces and envelopes

```python
from pyquefts import analysis as an, plot_surface, plot_envelope

biom = dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120)   # same as in q above
rs = an.response_surface(np.arange(0, 201, 20), np.arange(0, 61, 5), [60], "maize", soil=sup.to_soil(), biom=biom)
plot_surface(rs, x="N", y="P")               # filled contours at K = 60

env = an.uptake_yield_envelope("maize", biom=biom)
plot_envelope(env, res)                      # Janssen (1990) Fig. 1 lines, model runs overlaid
```

The envelope plot uses the yield of the modelled organ (`organ` output: leaves, stems or storage organ), so the overlaid points always lie between the maximum-accumulation and maximum-dilution lines.

## Inputs from a soil data cube

`run_points` and `run_spatial` give one result per location at a fixed fertilizer rate. For response curves the soil is extracted and Step 1 computed once, then the model is run for all rates:

```python
from pyquefts.geo import SoilCube, response_points, response_spatial, season_mean_temperature
from pyquefts import plot_curves, plot_spatial_response

cube = SoilCube.open("hnd_sh_soilgrids_format.nc", chunks={"depth": -1, "y": 1024, "x": 1024}, p_method="mehlich3")
temp = season_mean_temperature(start="05-15", end="08-31", bbox=(-89.4, 12.9, -83.1, 16.6))

# points: one curve per location (all points x all rates in one vectorized call)
df = response_points(cube, lon=[-87.65, -86.2], lat=[14.45, 14.8], crop="maize", nutrient="N",
                     rates=np.arange(0, 201, 20), temp=temp, store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                     fert=dict(N=0, P=0, K=0))
plot_curves(df)                                   # store_lim vs N, one line per "lon, lat" label
plot_curves(df, y=["store_lim", "N_gap"])
df.pivot(index="N", columns="label", values="store_lim")   # the table behind the plot

# region: run_spatial stacked along a rate dimension (dask-lazy until compute)
ds = response_spatial(cube, "maize", "N", [0, 40, 80, 120, 160], store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                      fert=dict(N=0, P=15, K=30), temp=temp, bbox=(-87.8, 14.3, -87.6, 14.5), compute=True)
plot_spatial_response(ds)                         # area mean with a 10-90 % quantile band
ds["store_lim"].plot(col="N", col_wrap=5)         # xarray facet maps, one per rate
```

`response_points` returns a long table (`point`, `label`, `lon`, `lat`, soil inputs, supplies, the rate column named after the nutrient, the two fixed rates as `P_applied`/`K_applied`, and the outputs). `response_spatial` returns the `run_spatial` Dataset with an extra dimension (`N`, `P` or `K`, kg/ha) and the supply maps; it can be written with `write_netcdf`. Both keep all `run_points`/`run_spatial` options (`temp`, `supply_method`, `bbox`, `mask`, depth layer, ...).

All functions return the matplotlib Axes and accept `ax=`/`axes=` to draw into an existing figure; extra keyword arguments go to `Axes.plot`.

`examples/04_response_plots.py` produces the five figures described above. Tests are in `tests/test_plot.py` (Agg backend; they check the data placed on the axes).
