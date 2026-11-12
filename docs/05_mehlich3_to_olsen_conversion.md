# Converting Mehlich-3 P to Olsen P for QUEFTS

## Why a conversion is needed

QUEFTS Step 1 (Janssen et al., 1990, Table II) relates the potential P supply to **P-Olsen** (0.5 M NaHCO3, pH 8.5; see `docs/01_scientific_review.md`, section 2) and the Sattari et al. (2014) recalibration keeps that method. The Suelos de Honduras soil survey, which provides the `p` variable of `hnd_sh_soilgrids_format.nc`, reports **Mehlich-3 extractable P** (confirmed by the data owner on 12 Nov 2025; the file metadata only says "mg/kg"). Mehlich-3 (acid fluoride extractant) typically extracts more P than Olsen on acid to neutral soils, so feeding the raw values to the model would overestimate P supply and underestimate P gaps. The values must be converted before Step 1.

R-QUEFTS has no such conversion; this is a documented **extension** of PyQUEFTS. The model equations themselves are unchanged: the conversion is applied to the input only.

## Equation implemented

From McDowell, Pletnyakov, Noble & Haygarth (2023), *A global database of soil plant available phosphorus*, Scientific Data 10:125, Table 2 and "Technical Validation" section. The authors did not fit new regressions; they computed sample-size-weighted means of published linear regressions (all in mg/kg):

| Source | Soils | pH range | n | Slope | Intercept | R² |
|---|---|---|---|---|---|---|
| Wolf & Baker (1985), USA | non-calcareous | 4.3-6.8 | 91 | 0.33 | 6.9 | 0.88 |
| Mallarino (1995), USA | non-calcareous | 5.3-8.2 | ~180 | 0.47 | 1.2 | 0.79 |
| Mallarino (1995), USA | calcareous | 5.3-8.2 | ~60 | 0.45 | 1.8 | 0.81 |
| Buondonno et al. (1992), Italy | non-calcareous | 3.8-8.6 | 66 | 0.70 | -0.6 | 0.71 |
| Buondonno et al. (1992), Italy | calcareous | 3.8-8.6 | 54 | 0.37 | 0.4 | 0.73 |
| **Weighted mean** | **non-calcareous** | | **337** | **0.47** | **2.4** | **0.80** |
| **Weighted mean** | **calcareous** | | **114** | **0.41** | **1.1** | **0.77** |

Implemented in `pyquefts/core/conversions.py`:

```
Olsen P = 0.47 * Mehlich-3 P + 2.4     non-calcareous soils (default)
Olsen P = 0.41 * Mehlich-3 P + 1.1     calcareous soils (calcareous=True)
```
Negative results are clamped to 0; NaN propagates. McDowell et al. validated the equations on 97 independent USDA-NCSS samples with paired measurements: regression slope of measured on predicted Olsen P 0.93, Nash-Sutcliffe efficiency above 0.7, and concluded that "the equations are suitable for general use as a conversion tool". They also note, citing Steinfurth et al. (2021), that pH has a minimal effect on Mehlich-3 to Olsen conversions, which is why separate equations are kept only for calcareous vs non-calcareous soils.

**Why the non-calcareous equation for Honduras.** The topsoil pH in the survey is mostly 5.5-7.3 (see the cube inspection in `docs/03_architecture.md`); calcareous horizons are the exception. The default is therefore `calcareous=False`; a per-pixel boolean (for example derived from a soil-class map) can be passed to `mehlich3_to_olsen` where a calcareous mask is available.

## Alternatives considered

* **Breure, Van Eynde, Kempen, Comans & Hoffland (2022)**, *Transfer functions for phosphorus and potassium soil tests and implications for the QUEFTS model*, Geoderma 406:115458 (doi 10.1016/j.geoderma.2021.115458). This is the only transfer function derived specifically for QUEFTS on tropical soils (Rwanda / sub-Saharan Africa). Its P function is ln(P-Olsen) = 0.769 ln(P-M3) + 0.620 ln(Al-M3): it requires **Mehlich-3 extractable Al**, which the Honduras cube does not contain, so it cannot be applied here. If Al-M3 becomes available it should be preferred and can be added to `P_METHODS` as a two-predictor method.
* **Steinfurth, Hirte, Morel & Buczko (2021)**, *Conversion equations between Olsen-P and other methods used to assess plant available soil phosphorus in Europe - A review*, Geoderma 401:115339. Reports Mehlich-3/Olsen ratios of about 1.6 on CaCO3-free acid soils vs about 2.5 on calcareous soils and stresses that coefficients from different studies differ by factors above two. Consistent with the 0.47 slope used (1/0.47 = 2.1 including the intercept effect at typical levels) but gives no single equation preferable to McDowell's weighted mean.
* **Dari et al. (2019)**, Idaho soils (SSSAJ 83:817), Olsen = Mehlich-3/2.05 - 14: fitted on calcareous and neutral soils of the US Pacific Northwest; the large negative intercept makes it unsuitable at the low P levels typical of Honduran soils (it returns 0 below about 29 mg/kg Mehlich-3).
* A 2026 EGUsphere preprint ("Unified conversion equations between Olsen and Mehlich-3 soil phosphorus tests") models the Olsen/M3 ratio as a function of soil properties. It is unreviewed and was not used.

## Effect on the Honduras cube

Typical Mehlich-3 P in the Comayagua test window is about 25 mg/kg (range 12-51). After conversion the Olsen estimate is about 14 mg/kg (8-26), inside the 0-30 mg/kg validity range of Janssen (1990). With `nut_supply2` the P supply also depends on total P, which is estimated as `ptotal_ratio * Olsen P` (default 55, Rquefts documentation; 95 per Sattari et al. 2014); the ratio is applied to the converted value.

## How to use

```python
from pyquefts.geo import SoilCube
cube = SoilCube.open("hnd_sh_soilgrids_format.nc", chunks={"depth": -1, "y": 1024, "x": 1024},
                     p_method="mehlich3")          # default is "olsen" (no conversion)
print(cube)                                        # shows "P method: mehlich3 -> Olsen (McDowell et al. 2023 ...)"

import pyquefts as pq
pq.run_point(ph=6.2, soc=15, kex=9, polsen=25, temp=24, p_method="mehlich3", crop="maize", ...)
pq.mehlich3_to_olsen(25)                           # 14.15
```
The output datasets of `supply_map`/`run_spatial` record the method in `attrs["p_method"]`.

## Uncertainty to keep in mind

The conversion carries the scatter of the underlying regressions (R² about 0.8) and was derived on temperate soils; no Central American calibration exists in the sources found. P supply, P uptake and the P gap inherit this uncertainty. A local Olsen vs Mehlich-3 paired dataset for Honduran soils would allow replacing the coefficients in `MEHLICH3_TO_OLSEN` without touching the rest of the model.

## References

* Buondonno, A., Coppola, E., Felleca, D., Violante, P., 1992. Comparing tests for soil fertility: 1. Conversion equations between Olsen and Mehlich 3 as phosphorus extractants for 120 soils of south Italy. Commun. Soil Sci. Plant Anal. 23, 699-716.
* Breure, T.S., Van Eynde, E., Kempen, B., Comans, R.N.J., Hoffland, E., 2022. Transfer functions for phosphorus and potassium soil tests and implications for the QUEFTS model. Geoderma 406, 115458.
* Dari, B., Rogers, C.W., Leytem, A.B., Schroeder, K.L., 2019. Evaluation of soil test phosphorus extractants in Idaho soils. Soil Sci. Soc. Am. J. 83, 817-824.
* Mallarino, A.P., 1995. Comparison of Mehlich-3, Olsen, and Bray-P1 procedures for phosphorus in calcareous soils. Proc. 25th North Central Extension-Industry Soil Fertility Conf., 96-101.
* McDowell, R.W., Pletnyakov, P., Noble, A., Haygarth, P.M., 2023. A global database of soil plant available phosphorus. Scientific Data 10, 125. https://doi.org/10.1038/s41597-023-02022-4
* Steinfurth, K., Hirte, J., Morel, C., Buczko, U., 2021. Conversion equations between Olsen-P and other methods used to assess plant available soil phosphorus in Europe - A review. Geoderma 401, 115339.
* Wolf, A.M., Baker, D.E., 1985. Comparisons of soil test phosphorus by Olsen, Bray I, Mehlich I and Mehlich III methods. Commun. Soil Sci. Plant Anal. 16, 467-484.
