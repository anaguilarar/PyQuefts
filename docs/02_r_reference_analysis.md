# Step 2 — Analysis of the reference implementation (R-QUEFTS / `Rquefts` 1.2-8)

Source inspected: GitHub `cropmodels/Rquefts` (identical to the CRAN 1.2-8 tarball, verified with a whitespace-insensitive diff). The model core is `src/quefts.cpp` (≈180 lines of C++); everything else is R glue.

## 1. Package layout

| File | Content | PyQUEFTS counterpart |
|---|---|---|
| `src/quefts.h` | structs `QueftsSoil`, `QueftsCrop`, `QueftsModel` | `pyquefts.core.parameters` dataclasses |
| `src/quefts.cpp` | `uptake()`, `yield()`, `requirements()`, `QueftsModel::run()`, `predict()`, `batch()` | `pyquefts.core.model` (vectorized NumPy) |
| `src/qUtil.h` | `approx()` piecewise-linear interpolation | `numpy.interp` with clamping |
| `R/nutSupply.R` | `nutSupply1`, `nutSupply2`, pH interpolation tables | `pyquefts.core.supply` |
| `R/rquefts.R` | defaults `quefts_soil()`, `quefts_crop()`, `quefts_fert()`, `quefts_biom()`; model constructor | `pyquefts.core.parameters`, `pyquefts.crops` |
| `R/ferti_funs.R` | `fertilizers()`, `nutrientRates()`, `fertApp()` (LP, limSolve), `optApp()` (L-BFGS-B on net revenue) | `pyquefts.fertilizer` (scipy `linprog`, `minimize`) |
| `R/reverse.R` | `revSupply()` — fit base supplies to omission-trial yields | `pyquefts.reverse` |
| `R/spatial.R` | `predict()` on terra SpatRaster, `batch()` on matrices | `pyquefts.geo.spatial`, `pyquefts.point.run_batch` |
| `inst/extdata/quefts_crop_pars.csv` | 28 crops × 14 parameters | copied verbatim to `pyquefts/data/` |
| `inst/extdata/fertilizers.csv` | 29 products, nutrient % | copied verbatim |

## 2. Data structures

### Soil (`QueftsSoil`)
```
N_base_supply, P_base_supply, K_base_supply   kg/ha, 120-day season      defaults 60, 10, 60
N_recovery, P_recovery, K_recovery             fraction                    defaults 0.5, 0.1, 0.5
UptakeAdjust                                   2-col matrix (days, factor) default (0,0)(40,.4)(80,.7)(120,1)(240,1.6)(360,2)(1000,2)
```
### Crop (`QueftsCrop`)
```
NminStore NminVeg NmaxStore NmaxVeg    kg/kg
PminStore PminVeg PmaxStore PmaxVeg
KminStore KminVeg KmaxStore KmaxVeg
Yzero                                  kg/ha vegetative biomass at zero storage yield
Nfix                                   fraction of N from fixation
```
Default (no name) = maize values. Named crops from the CSV (case-insensitive lookup).

### Biomass / management (fields of `QueftsModel`)
```
leaf_att, stem_att, store_att   attainable DM (kg/ha)      defaults 2500, 2500, 5000
SeasonLength                    days                        default 120
N, P, K                         fertilizer applied (kg/ha)  defaults 60, 10, 60
```

## 3. Calculation sequence in `QueftsModel::run()`

1. **Organ selection.** If `store_att > 200` → model the storage organ (`organs = 3`, `Yo = store_att`) and compute `Nzero, Pzero, Kzero` from `Yzero`. Else if `stem_att <= 200` → leaves (`organs = 1`, `Yo = max(5, leaf_att)`); else stems (`organs = 2`, `Yo = stem_att`). In the last two cases the zero-uptakes are 0.
2. **Season factor.** `relgrowseason = approx(UptakeAdjust, SeasonLength)`.
3. **Supply.** `X_supply = X_base_supply * relgrowseason + X_fertilizer * X_recovery`.
4. **Requirements.** `requirements()` returns `(reqAcc, reqDil, reqFert)` per nutrient; `X_gap = reqFert`.
5. **Yield/uptake ratios.** `XYratA = Yo / max(0.1, reqAcc - Xzero)`, `XYratD = Yo / max(0.1, reqDil - Xzero)`.
6. **Uptake.** Six pairwise `uptake()` calls (Table IV with the crop ratios); `UX = max(0, min(pair1, pair2))`.
7. **Yield ranges.** `YXA = XYratA * max(0, UX - Xzero)`, `YXD = XYratD * max(0, UX - Xzero)`.
8. **Pairwise yields.** Six `yield()` calls (Table V with boundary rules) → EYNP, EYPN, EYNK, EYKN, EYPK, EYKP.
9. **Final yield.** `min(min(YND, YPD, YKD, Yo), mean(six))`.
10. **Organ partition.** For `organs == 3`: `store_lim = yield`; `leaf_lim = leaf_att/(leaf_att+stem_att) * ((stem_att+leaf_att - Yzero) * yield/store_att + Yzero)`; `stem_lim = stem_att/leaf_att * leaf_lim`; with a fallback when `yield <= 0`. For organs 1 and 2 the other organs scale proportionally.

## 4. Control flow of `yield()` (needed for exact vectorization)

```
LYD = min(iYD, jYD, zYD, Yo)
if iYD <= jYA:                 -> return LYD           (iYD > LYD is impossible, so EY == LYD)
elif iYA >= LYD or iYD <= LYD:
    if iYA < LYD and iYD == LYD:  LYD = min(jYD, zYD, Yo)   # fall through to the parabola
    else:                         -> return LYD
EY = parabola(jYA, LYD, iUPT - izero, iYratA, iYratD)
if jYA > LYD: EY = LYD
return EY
```
Since `LYD <= iYD` always holds, the condition `iYD <= LYD` is equivalent to `iYD == LYD`. The branch structure reduces to three masks that are implemented with `numpy.where` in PyQUEFTS.

## 5. Interfaces reproduced

| R call | PyQUEFTS |
|---|---|
| `nutSupply1(pH, SOC, Kex, Polsen, Ptotal=NA)` | `nut_supply1(ph, soc, kex, polsen, ptotal=None)` |
| `nutSupply2(temp, pH, SOC, Kex, Polsen, Ptotal)` | `nut_supply2(temp, ph, soc, kex, polsen, ptotal)` |
| `quefts(soil, crop, fert, biom); run(q)` | `Quefts(soil, crop, fert, biom).run()` → `QueftsResult` |
| `batch(q, supply, fert, yatt, leaf_ratio, stem_ratio, var)` | `run_batch(...)` (vectorized) |
| `predict(q, supply_raster, yatt_raster, ...)` | `pyquefts.geo.spatial.run_spatial(...)` (xarray/dask) |
| `fertilizers(); nutrientRates(); fertApp(); optApp()` | `pyquefts.fertilizer` |
| `revSupply(obs, crop, soil, Ya, leaf_ratio, stem_ratio)` | `pyquefts.reverse.rev_supply(...)` |
| `quefts_crop(name)` | `pyquefts.crops.get_crop(name)` |

## 6. Reference outputs captured for validation

Documented in the reagro.org tutorial and the package examples (all reproduced by the installed R package 1.2-8 on this machine):

* default model `run(quefts())`: store_lim 3756.94929, N_uptake 81.66872, P_uptake 10.88951, K_uptake 77.38765, N_gap 60, P_gap 205, K_gap 0.
* `nutSupply2(25, 6, 15, 1, 10, 500)` → N 76.03381, P 12, K 18.32258.
* Barley, N=100, biomass (2200, 2700, 4800), season 110 → store_lim 2877.4, N_supply 105.5, P_gap 160.9, K_gap 113.3.
* Rice N-response series (leaf 2651, stem 5053, store 8208): 3162.406 at N=1 … 4384.379 at N=191.

The validation framework (`validation/`) additionally generates several thousand randomized R runs and compares them with PyQUEFTS at double precision.
