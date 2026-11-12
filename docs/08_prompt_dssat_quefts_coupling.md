# Prompt: couple WeatherSoilDataProcessor (DSSAT) outputs with PyQUEFTS for probabilistic nutrient recommendations

Paste everything below the line into the agent that works in the WeatherSoilDataProcessor project.

---

## Role and goal

You are extending the project `D:\OneDrive - CGIAR\scripts\SueloDeHonduras\WeatherSoilDataProcessor` (village-level crop simulation for Honduras with DSSAT) with a coupling module that turns its multi-year, multi-planting-date yield simulations into **probabilistic QUEFTS nutrient-limited yields and N, P, K recommendations**. QUEFTS is provided by the installed Python package **PyQUEFTS** (`D:\OneDrive - CGIAR\scripts\SueloDeHonduras\PyQuefts`, importable as `pyquefts` in the conda env `gs_env`). You must not re-implement or modify QUEFTS equations; you call PyQUEFTS. The user makes all git commits; never commit.

The scientific logic of the coupling, which you must follow:

1. QUEFTS is a static, deterministic nutrient balance. Weather enters it only through (a) the **attainable yield and biomass of the crop** (`store_att`, `leaf_att`, `stem_att`), (b) the **season length** (days, through the crop's uptake-adjust table) and (c) the **growing-season mean temperature** (Sattari et al. 2014 N supply). DSSAT provides all three per year and planting date.
2. Therefore the probabilistic version is obtained by running the unchanged QUEFTS once per (soil pixel x year x planting date) and summarizing the resulting distribution, not by changing the model.
3. Soil chemistry (pH, SOC, Olsen P, exchangeable K) varies pixel by pixel inside a texture class. Do not average soil inputs per texture: the model is non-linear (quadratic pH factors, parabolic uptake-yield curves), so the output of the mean is not the mean of the outputs. Keep pixels and attach to each the yield ensemble of its texture class.

## What exists (verified by inspection; re-verify before relying on it)

Run folder example: `runs/070331/` (village GEOCODIGO 070331, bean).

* `bean_potential_yield.csv`: 2800 rows = 35 years (1990-2024) x 40 weekly planting dates (TRNO 1..40, 1 Feb to 31 Oct) x 2 texture classes (`loam`, `sandyclayloam`). Columns: index, `PDAT`, `MDAT`, `HDAT` (dates), `HWAH` (harvested yield, kg/ha **dry matter**, DSSAT convention), `TRNO`, `texture`. Yields range 2-3039 kg/ha; February-March plantings give ~100-300 kg/ha, May-August ~2200 kg/ha.
* `loam/Summary_<TRNO>.OUT`, `sandyclayloam/Summary_<TRNO>.OUT` (40 each): DSSAT summary with `CWAM` (tops biomass at maturity, kg/ha DM), `HWAM`, `HWAH`, `BWAH`, `HIAM`, `PDAT`, `ADAT`, `MDAT`, `HDAT`, `PRCM` (season rainfall), `NFXM` (N fixed), etc. Use these to get **vegetative biomass = CWAM - HWAM** per run; `bean_potential_yield.csv` alone lacks biomass.
* `loam/EXPS0001.BNX` simulation controls: `WATER Y, NITRO Y, SYMBI Y, PHOSP N, POTAS N`; 35 years; fertilizer 200 kg N at 30 days after planting (`model_configurations/crop_model_configuration_bean.yaml`). **Consequence:** these yields are water-limited and N-simulated under heavy N; they are not strictly "N-unlimited". Feeding them to QUEFTS as attainable yield is acceptable only if N was effectively non-limiting in DSSAT. Check this (compare with a run at `NITRO N`, or inspect N stress in PlantGro/Summary); if N stress is material, re-run DSSAT with `NITRO N` (water-limited potential) for the attainable yields. Document the choice. P and K are not simulated by DSSAT, so there is no double counting for them.
* `weatherbytreatment.csv`: daily `TMIN`, `TMAX`, `RAIN` per `TRNO` and `DATE` (the DSSAT weather actually used). Season mean temperature per run = mean of (TMIN+TMAX)/2 from PDAT to MDAT.
* `soil_.nc`: village soil cube on a 0.005 deg grid (17 x 15 cells, EPSG:4326, nodata -99.99, `scale_values_factor: 100` attribute, check what it means), depths `0-5, 5-15, 15-30, 30-60`, variables `sand, clay, silt, bdod, phh2o, soc, cec, k, p, wv0010, wv0033, wv1500, cfvo, nitrogen`. Observed top-layer means: phh2o ~59 (pH x 10), soc ~200 (dg/kg), p ~21 (mg/kg, **Mehlich-3**: the Honduras survey P; confirm with the user), nitrogen ~200 (cg/kg), k (expected cmol(+)/kg, verify). These follow the Honduras cube conventions handled by `pyquefts.geo.SoilCube` (`UNIT_FACTORS`), but `soil_.nc` may lack `units` attributes, so pass an explicit mapping.
* `texture.tif`: USDA texture class codes per cell (here 4 and 7 for 9 and 171 cells). Find the code-to-name mapping in `spatialdata/soil_data.py` (the same one that produced `loam` / `sandyclayloam`), do not guess it.
* `soil_and_nitrogen_values_bean.json`: per-texture mean phh2o, SOM (%) and total N used to build the DSSAT soil profile. Do not use it for QUEFTS inputs (means, and SOM not SOC); use the pixel values of `soil_.nc`.
* Village polygons: `data/tb_limitealdeas.shp`, key `GEOCODIGO`; crop areas per village in `data/*_area_village_spam2020.csv`.

PyQUEFTS facts you need:

* Crop table: `pyquefts.get_crop("Field bean")` (Nfix = 0.75), `get_crop("maize")`, `list_crops()`. **Quirk inherited from R-QUEFTS, documented in `PyQuefts/docs/01_scientific_review.md` section 10:** when a crop has `Nfix != 0` the model reports `N_gap = 0` always. For beans the N recommendation from QUEFTS is therefore not meaningful; report P and K gaps for beans and say so explicitly. For maize all three are valid.
* Organ selection: `store_att > 200` kg/ha means the storage organ is modelled. Runs with HWAH <= 200 (February plantings, failed seasons) switch QUEFTS to modelling stems or leaves. Keep them in the ensemble (they are real failure years) but flag them with the `organ` output; do not average them with grain yields without saying so.
* Biomass: `pyquefts.Biomass(leaf_att, stem_att, store_att, SeasonLength)` or `Biomass.from_yield(store_att, leaf_ratio, stem_ratio, SeasonLength=...)`. Prefer DSSAT biomass: `store_att = HWAH`, `leaf_att + stem_att = CWAM - HWAM`, split leaf:stem with a documented ratio (DSSAT `PlantGro.OUT` columns `LWAD`/`SWAD` at maturity if output is enabled, otherwise a cited literature ratio for common bean; state the source). `SeasonLength = MDAT - PDAT` in days (physiological maturity; justify if you use HDAT).
* Soil supply: `pyquefts.nut_supply2(temp, ph, soc, kex, polsen, ptotal=None, ptotal_ratio=None)` (Sattari 2014; array inputs broadcast); `pyquefts.geo.SoilCube(ds, mapping=..., p_method="mehlich3")`, `cube.model_inputs(["ph","soc","kex","polsen"], top=0, bottom=20)` gives thickness-weighted 0-20 cm inputs in model units; `pyquefts.geo.supply_map(cube, temp=...)`. Mehlich-3 P is converted to Olsen P with McDowell et al. 2023 (`PyQuefts/docs/05_mehlich3_to_olsen_conversion.md`).
* Model: `pyquefts.run_quefts(soil, crop, fert, biom)` and `pyquefts.run_batch(crop, supply(n,3), fert(n,3), yatt, leaf_ratio, stem_ratio, soil=, season_length=)` are fully vectorized over NumPy arrays of any shape: one call for all pixels x years x planting dates is the intended use (the core processes ~300k-500k elements per second). Outputs: `leaf_lim, stem_lim, store_lim, N/P/K_supply, N/P/K_uptake, N/P/K_gap` plus intermediates (`organ`, `YND`, ...). `pyquefts.analysis.limiting_nutrient(res)`.
* Fertilizer economics: `pyquefts.fertilizer.opt_app`, `fert_app`, `nutrient_rates`, `fertilizers()`; Honduras prices may be in `data/fertilizer_list_hnd.csv`.
* Plots: `pyquefts.plot` (`plot_result`, `plot_curves`, `plot_response_curves`, `plot_spatial_response`), see `PyQuefts/docs/07_plotting.md`.
* Read `PyQuefts/README.md`, `docs/03_architecture.md` and `docs/07_plotting.md` before writing code.

## Tasks

Work step by step; after each step report what you verified, with numbers.

**Step 0, verification.** Confirm the DSSAT switches and fertilizer, HWAH units, the texture code mapping, the units and nodata of `soil_.nc`, that `p` is Mehlich-3, and that `k` is exchangeable K in cmol(+)/kg. Decide (with the user if needed) whether the current `NITRO Y` yields are acceptable as attainable yields or a `NITRO N` re-run is required. Do not proceed on unverified assumptions.

**Step 1, attainable-yield ensemble.** Parse the 80 `Summary_*.OUT` files into one table per village with columns `texture, TRNO, year, PDAT, MDAT, HDAT, season_length, HWAH, CWAM, veg_biomass, leaf_att, stem_att, store_att, temp_season, rain_season`. Join the season temperature from `weatherbytreatment.csv`. Check row count = textures x 40 x 35 and that HWAH matches `bean_potential_yield.csv`.

**Step 2, soil inputs per pixel.** Open `soil_.nc` with `SoilCube` and an explicit mapping (units factors as above, nodata handled, `p_method="mehlich3"`), aggregate 0-20 cm, and produce a per-pixel table `pixel_id, x, y, texture_class, ph, soc, kex, polsen`. Drop nodata pixels and report how many remain per texture.

**Step 3, probabilistic QUEFTS run.** For each texture class, form the cross product pixels x (year, TRNO) and run QUEFTS in one vectorized call per fertilizer scenario: supply from `nut_supply2` with the run's season temperature, biomass from Step 1, season length from Step 1. Scenarios: no fertilizer, and a grid of rates (for beans P and K only, e.g. P 0-60 kg/ha, K 0-120 kg/ha; for maize N 0-200, P 0-60, K 0-120). Store the full result as a tidy table (parquet) and as an xarray Dataset with dims (`pixel`, `year`, `TRNO`, scenario dims).

**Step 4, summaries that make the recommendation probabilistic.** Per village, texture and planting date (TRNO) and for the village as a whole (area-weighted by pixel counts, optionally by crop area from SPAM):

* distribution of nutrient-limited yield `store_lim` (mean, p10, p25, p50, p75, p90) versus the attainable `store_att`, i.e. the nutrient yield gap distribution;
* probability that `store_lim` exceeds thresholds (e.g. 50 %, 75 %, 90 % of attainable);
* distribution of P and K gaps (and N for maize) and of the limiting nutrient;
* **risk-based recommendation**: for each nutrient the lowest rate at which the model yield reaches at least X % of the attainable yield in at least Y % of the years (defaults X = 90, Y = 75; make both parameters), and, if prices are available, the rate maximizing expected profit over the ensemble via `opt_app` per year (or a grid), reporting also the profit quantiles;
* a comparison of these with the deterministic recommendation obtained from the 35-year mean attainable yield, to show the difference the ensemble makes.

**Step 5, optional uncertainty layer (only after Steps 1-4 work).** Monte Carlo over input uncertainty, each source documented with its reference: residual scatter of the Mehlich-3 to Olsen regression (McDowell et al. 2023, R2 0.80), the total-P ratio (`ptotal_ratio` 55 vs 95, see `PyQuefts/docs/01_scientific_review.md`), and literature ranges of the fertilizer recoveries (`SoilParameters(N_recovery, P_recovery, K_recovery)`). Report a simple variance decomposition (which source dominates the yield and the recommendation).

**Step 6, outputs and plots.** Write per-village CSV/parquet summaries and a NetCDF of pixel-level quantiles; produce figures: yield distributions by planting date (box or fan plots) for attainable vs nutrient-limited; P and K response curves with year-to-year bands (`plot_curves` / `plot_spatial_response` style); maps of the p50 nutrient-limited yield and of the recommended rates; the deterministic vs probabilistic comparison. Write a short methods document stating every assumption, unit, equation source and the Nfix limitation for beans.

**Step 7, tests.** Add pytest tests: the ensemble run for one pixel and one (year, TRNO) equals a direct `pyquefts.Quefts(...).run()` with the same inputs; row counts; nodata handling; the Jensen check (mean of per-pixel outputs differs from output of the mean inputs); recommendation monotonicity (higher Y requires a rate at least as high).

## Rules

* Do not change QUEFTS equations or parameters; everything probabilistic comes from inputs and post-processing. If you believe a model change is needed, stop and explain instead.
* Do not invent units, conversions, ratios or parameter values; if a value is not in the data or a cited source, ask.
* Keep DSSAT dry-matter units throughout (QUEFTS also uses kg/ha dry matter). If farmer-moisture yields are needed for reporting, convert at the end and say so.
* Vectorize with NumPy/xarray; no Python loops over pixels or years.
* Report faithfully: if a check fails or a step is skipped, say so with the numbers.
