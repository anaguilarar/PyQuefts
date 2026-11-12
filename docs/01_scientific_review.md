# Step 1 — Scientific review of QUEFTS

Sources reviewed (project `bibliography/` folder):

| Reference | Role in this project |
|---|---|
| Janssen et al. (1990), *Geoderma* 46:299–318 | Original model. Defines the four calculation steps, Table II (soil → supply), Tables III/IV (supply → uptake), Table V (uptake → yield, pairwise combination). |
| Smaling & Janssen (1993), *Geoderma* 59:21–44 ("Calibration of QUEFTS") | Recalibration for Kenya; introduces fertilizer recovery fractions, a harvest-index constraint and a modified uptake curve. R-QUEFTS keeps the 1990 uptake/yield equations. |
| Sattari et al. (2014), *Field Crops Research* 157:35–46 | Tests QUEFTS in several environments. Table 4 gives recalibrated Step-1 parameters (temperature effect on N mineralization, piecewise pH functions, Ptotal ≈ 95·P-Olsen). Basis of `nutSupply2()` in R-QUEFTS. |
| Ravensbergen et al. (2021), *FCR* 266:108126 | Quantile-regression calibration of Step 1 for smallholder data. Not implemented; relevant for future calibration work. |
| Yang et al. (2017), *PLoS ONE* 12:e0177509 | QUEFTS used to derive balanced nutrient-uptake requirements for soybean (the "internal efficiency" use of the model). Supports the balanced-requirement output. |
| Rquefts manual (CRAN v1.2-8) and reagro.org tutorial | Reference implementation and documented example outputs. |

## 1. What QUEFTS computes

QUEFTS treats soil fertility as the capacity of the 0–20 cm topsoil to supply N, P and K and expresses it as the yield of a reference crop. It is a static, season-scale model with four successive steps (Janssen 1990, "Framework of the system"):

1. **Potential supply** of N, P and K (kg/ha) from a few routine topsoil properties.
2. **Actual uptake** of each nutrient as a function of its own potential supply and the supplies of the other two (law of the minimum, applied pairwise).
3. **Yield ranges** for each nutrient from its uptake: from maximum accumulation (lower bound) to maximum dilution (upper bound).
4. **Ultimate yield**: combine the yield ranges in pairs (six parabolas), cap each by the dilution limit of the third nutrient and by the attainable yield, and average the six estimates.

Validity domain of Step 1 (Janssen 1990): deep, well-drained soils; pH(H2O) 4.5–7.0; organic C < 70 g/kg; P-Olsen < 30 mg/kg; exchangeable K < 30 mmol/kg; no other limiting growth factors. Steps 2–4 are theoretically derived and considered more generally applicable.

## 2. Step 1 — Soil properties to potential supply

### 2.1 Original equations (Janssen 1990, Table II; R-QUEFTS `nutSupply1`)

Units: organic C (SOC) in g/kg; P-Olsen and total P in mg/kg; exchangeable K in mmol/kg; pH in water; sample depth 0–20 cm. Output in kg/ha per standard growing season.

```
fN = 0.25 (pH - 3)                        (eq. 4)
fP = 1 - 0.5 (pH - 6)^2                   (eq. 5)
fK = 0.625 (3.4 - 0.4 pH)                 (eq. 6)

SN = fN * 6.8 * SOC                       (eq. 1)   [or fN * 68 * Norg, eq. 7, assumes C:N = 10]
SP = fP * 0.35 * SOC + 0.5 * Polsen       (eq. 2)
SP = fP * 0.014 * Ptotal + 0.5 * Polsen   (eq. 8, preferred when total P is known)
SK = fK * 400 * Kex / (2 + 0.9 * SOC)     (eq. 3)
```

Eq. 2 assumes Ptotal/SOC ≈ 25 (mg/kg per g/kg), which holds in unfertilized soils. R-QUEFTS `nutSupply1` uses eq. 2 unless `Ptotal` is given, and clamps negative supplies to 0.

### 2.2 Recalibrated equations (Sattari et al. 2014, Table 4; R-QUEFTS `nutSupply2`)

* The N mineralization coefficient becomes temperature dependent: `alphaN = 2 * 2^((T - 9)/9)` replaces 6.8 (T = mean growing-season temperature in °C; equals 6.8 at T ≈ 24.9 °C).
* P supply always uses total P: `SP = fP * 0.014 * Ptotal + 0.5 * Polsen`. When total P is not measured, Sattari et al. suggest `Ptotal ≈ 95 * Polsen` for unfertilized soils; the Rquefts author reports `≈ 55 * Polsen` from AfSIS data; Janssen's `25 * SOC` is a third option.
* The pH factors become piecewise. R-QUEFTS implements them as piecewise-linear interpolation tables (exact for fN, close approximations for fP and fK):

```
fN: (4.7, 0.4) (7, 1)                                                        constant outside
fP: (4.7, 0.02) (5, 0.5) (5.5, 0.875) (6, 1) (6.7, 1) (7, 0.9775) (7.5, 0.84) (8, 0.02)   constant outside
fK: (4.5, 1) (5, 0.8842312) (6, 0.71) (6.8, 0.6)                             constant outside
```

The K equation keeps alphaK = 400 (Sattari's recalibrated 500 is **not** used by R-QUEFTS). `nutSupply2` does not clamp negatives.

### 2.3 Season length

Supplies refer to a standard 120-day season. R-QUEFTS rescales the base supply with a piecewise-linear `UptakeAdjust` table (days → multiplier): (0,0) (40,0.4) (80,0.7) (120,1) (240,1.6) (360,2) (1000,2). The table is an R-QUEFTS addition, not part of Janssen 1990.

### 2.4 Fertilizer

Total supply = base supply × season factor + fertilizer × recovery fraction (Smaling & Janssen 1993; Sattari 2014 eqs. 1–3). Default recoveries RN = 0.5, RP = 0.1, RK = 0.5.

## 3. Step 2 — Potential supply to actual uptake (Janssen 1990, Table IV)

For nutrient 1 given nutrient 2, with yield/uptake slopes at accumulation `a` and dilution `d` and minimum uptake `r` (generalized in R-QUEFTS to crop-specific ratios, see section 5):

```
Situation A:  S1 <= r1 + (S2 - r2) * a2/d1                      ->  U1(2) = S1
Situation C:  S1 >= r1 + (S2 - r2) * (2*d2/a1 - a2/d1)          ->  U1(2) = r1 + (S2 - r2) * d2/a1
Situation B:  otherwise
              U1(2) = S1 - 0.25 * [S1 - r1 - (S2 - r2)*a2/d1]^2 / [(S2 - r2) * (d2/a1 - a2/d1)]
```

Actual uptake U1 = min(U1(2), U1(3)); clamped at 0 in R-QUEFTS.

## 4. Step 3 — Uptake to yield range (Janssen 1990, Table III)

```
Y1A = a1 * (U1 - r1)     maximum accumulation (lower bound)
Y1D = d1 * (U1 - r1)     maximum dilution (upper bound)
```
Maize constants: N a=30 d=70 r=5; P a=200 d=600 r=0.4; K a=30 d=120 r=2 (grain at 12 % moisture).

## 5. Generalization to crops and organs (R-QUEFTS; Janssen 1993 "QUEFTS for other crops")

R-QUEFTS replaces the maize constants with crop-specific minimum and maximum nutrient mass fractions in vegetative (`minVeg`, `maxVeg`) and storage organs (`minStore`, `maxStore`), plus the attainable dry-matter biomass of leaves, stems and storage organ. From these:

```
reqAcc = (maxVeg*(leaf+stem) + maxStore*store) * (1 - Nfix)   uptake at maximum accumulation
reqDil = (minVeg*(leaf+stem) + minStore*store) * (1 - Nfix)   uptake at maximum dilution
Xzero  = Yzero * (maxVeg + 2*minVeg)/3  [* (1 - Nfix) for N]  uptake at zero storage yield (the role of r)
YratA  = Yo / max(0.1, reqAcc - Xzero)                        the role of a
YratD  = Yo / max(0.1, reqDil - Xzero)                        the role of d
```
`Yo` is the attainable yield of the organ being modelled (storage organ when store_att > 200 kg/ha; otherwise stems, or leaves). `Nfix` is the fraction of N uptake from biological fixation (0.75 for legumes in the crop table). `Yzero` is the maximum vegetative biomass at zero storage-organ yield.

## 6. Step 4 — Combining yield ranges (Janssen 1990, Table V)

For nutrient pair (1, 2) the yield Y12 is a parabola in U1 that starts at Y2A (where Y1D = Y2A) and peaks at the upper limit L:

```
Y12 = Y2A + 2(L - Y2A)(U1 - r1 - Y2A/d1)/(L/a1 - Y2A/d1) - (L - Y2A)(U1 - r1 - Y2A/d1)^2/(L/a1 - Y2A/d1)^2
L   = min(Y1D, Y2D, Y3D, Yattainable)      Janssen: "if Y3D < Y2D substitute Y3D"
```
plus boundary rules: if Y1D <= Y2A the yield is the lower limit L; if Y1A >= L the yield is L. The final yield is min(Y1D, Y2D, Y3D, Yatt, mean of the six pairwise yields).

Worked example in Janssen 1990 Table VI (Soil A: pH 6, SOC 20, Polsen 5, Kex 10): SN 102, SP 9.5, SK 125 → UN 94, UP 9.3, UK 108 → YE 4279 kg/ha. The example uses the maize constants of Table III; the PyQUEFTS test suite reproduces it with an equivalent crop parameterization.

## 7. Fertilizer requirement ("gap")

R-QUEFTS computes, per nutrient, the fertilizer needed to reach the attainable biomass at *balanced* (mid-range) concentrations:
```
req = ((minVeg+maxVeg)/2 * (leaf+stem) + (minStore+maxStore)/2 * store)     (0 for N when Nfix > 0)
gap = max(0, (req - supply) / recovery)
```
This is the "balanced nutrient requirement" of Janssen (1998) and Yang et al. (2017) applied per nutrient. It ignores interactions, so it is an upper-bound recommendation for each nutrient in isolation.

## 8. Outputs required of PyQUEFTS (R-QUEFTS `run()` returns 12 values)

| Output | Meaning |
|---|---|
| `N_supply, P_supply, K_supply` | total supply (soil × season + fertilizer × recovery), kg/ha |
| `N_uptake, P_uptake, K_uptake` | actual uptake, kg/ha |
| `leaf_lim, stem_lim, store_lim` | nutrient-limited dry-matter yield of leaves, stems and storage organ, kg/ha |
| `N_gap, P_gap, K_gap` | fertilizer needed to reach the attainable biomass, kg/ha |

PyQUEFTS additionally exposes the intermediate quantities of each step (pairwise uptakes, yield ranges YxA/YxD, the six pairwise yields), because they are the basis of the nutrient–yield relationship outputs requested in the project brief.

## 9. Soil-test methods other than P-Olsen (extension)

Step 1 is calibrated on P-Olsen (0.5 M NaHCO3) and NH4-acetate exchangeable K. When a survey reports Mehlich-3 P, PyQUEFTS converts it to Olsen P before Step 1 using the sample-size-weighted regressions compiled by McDowell et al. (2023, *Scientific Data* 10:125, Table 2): Olsen = 0.47·M3 + 2.4 for non-calcareous soils (n = 337, R² 0.80) and 0.41·M3 + 1.1 for calcareous soils (n = 114, R² 0.77), validated on 97 USDA-NCSS samples (slope 0.93). Breure et al. (2022, *Geoderma* 406:115458) derived a QUEFTS-specific transfer function for tropical soils, ln(P-Olsen) = 0.769·ln(P-M3) + 0.620·ln(Al-M3), which needs Mehlich-3 Al and is therefore not applicable to soil maps without it. Converted values carry the regression uncertainty (coefficients across studies differ by factors of up to two, Steinfurth et al. 2021), which should be remembered when interpreting P supply and P gaps.

## 10. Known quirks in the reference implementation (documented; reproduced where harmless)

1. `approx()` in `qUtil.h` drops the last x-breakpoint of `UptakeAdjust`. Harmless for the default table (the last two y values are equal) but wrong for custom tables. PyQUEFTS uses a correct clamped linear interpolation; the validation suite confirms identical results for the default table.
2. When the computed yield is <= 0, the leaf-biomass fallback uses `(NminVeg + 2*NminVeg)` where `(NmaxVeg + 2*NminVeg)` was evidently intended (it would then match the `Nzero` definition). Reproduced for fidelity and flagged in the code.
3. In `requirements()` the N gap is forced to 0 for any crop with `Nfix > 0`, regardless of the fixation fraction.
4. `nutSupply1` clamps negative supplies to zero; `nutSupply2` does not.
