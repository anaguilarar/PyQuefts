# PyQUEFTS vs R-QUEFTS validation report

Reference: R version 4.4.1 (2024-06-14 ucrt); Rquefts 1.2.8; generated 2026-11-12 09:47:58
Tolerance: rtol=1e-09, atol=1e-07

## Model runs (5000 randomized scenarios, 29 crop parameter sets)

| variable   |    n |   n_mismatch |   max_abs_diff |   max_rel_diff |   n_nan_r |   n_nan_py |
|:-----------|-----:|-------------:|---------------:|---------------:|----------:|-----------:|
| leaf_lim   | 5000 |            0 |    4.91127e-11 |    4.62798e-15 |         0 |          0 |
| stem_lim   | 5000 |            0 |    4.91127e-11 |    4.87277e-15 |         0 |          0 |
| store_lim  | 5000 |            0 |    4.91127e-11 |    4.90782e-15 |         0 |          0 |
| N_supply   | 5000 |            0 |    5.00222e-12 |    4.88307e-15 |         0 |          0 |
| P_supply   | 5000 |            0 |    5.11591e-13 |    2.1199e-13  |         0 |          0 |
| K_supply   | 5000 |            0 |    5.00222e-12 |    4.70227e-15 |         0 |          0 |
| N_uptake   | 5000 |            0 |    4.54747e-12 |    4.89033e-15 |         0 |          0 |
| P_uptake   | 5000 |            0 |    4.83169e-13 |    2.1199e-13  |         0 |          0 |
| K_uptake   | 5000 |            0 |    4.54747e-12 |    1.58387e-14 |         0 |          0 |
| N_gap      | 5000 |            0 |    5.00222e-12 |    4.79102e-15 |         0 |          0 |
| P_gap      | 5000 |            0 |    4.77485e-12 |    4.91884e-15 |         0 |          0 |
| K_gap      | 5000 |            0 |    5.00222e-12 |    4.76247e-15 |         0 |          0 |

Scenarios with any mismatch: **0**

Scenario design: random base supplies (incl. zeros and extremes), recoveries, fertilizer rates (30 % zeros), attainable yields 500-15000 kg/ha with 8 % tiny storage organs (exercising the leaf/stem organ branches), season lengths 60-400 days, all 28 crops plus the default set.

## Supply functions (2000 random soils, pH 3.5-9.5, SOC 0-80, Kex 0-40, P-Olsen 0-40, Ptotal 0-2000, T 10-35)

| variable            |    n |   n_mismatch |   max_abs_diff |   max_rel_diff |
|:--------------------|-----:|-------------:|---------------:|---------------:|
| s1_N_base_supply    | 2000 |            0 |    1.59162e-12 |    1.11439e-14 |
| s1_P_base_supply    | 2000 |            0 |    2.46914e-13 |    4.57934e-12 |
| s1_K_base_supply    | 2000 |            0 |    1.00044e-11 |    4.73407e-13 |
| s1tot_N_base_supply | 2000 |            0 |    1.59162e-12 |    1.11439e-14 |
| s1tot_P_base_supply | 2000 |            0 |    2.23599e-13 |    4.35957e-13 |
| s1tot_K_base_supply | 2000 |            0 |    1.00044e-11 |    4.73407e-13 |
| s2_N_base_supply    | 2000 |            0 |    3.41061e-12 |    8.69918e-15 |
| s2_P_base_supply    | 2000 |            0 |    2.52243e-13 |    1.28115e-13 |
| s2_K_base_supply    | 2000 |            0 |    1.04592e-11 |    9.68817e-15 |

## Documented examples (reagro.org tutorial / Rquefts manual)

| case        |   max_abs_diff |
|:------------|---------------:|
| barley_N100 |    4.77485e-12 |
| barley_long |    1.81899e-12 |
| default     |    3.63798e-12 |
| rice_N1     |    4.09273e-12 |
| rice_N101   |    3.18323e-12 |
| rice_N191   |    1.81899e-12 |
| rice_N51    |    5.00222e-12 |
