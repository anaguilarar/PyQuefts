"""Soil-test conversions to the methods QUEFTS was calibrated on.

QUEFTS Step 1 (Janssen et al., 1990) uses P-Olsen (0.5 M NaHCO3) and
NH4-acetate exchangeable K. Soil surveys often report other extractants;
these helpers convert them with *published* regressions. They are an
extension to R-QUEFTS (which has no conversions) and every coefficient is
referenced so the choice can be revisited.

Phosphorus, Mehlich-3 -> Olsen (mg/kg)
--------------------------------------
McDowell, R.W., Pletnyakov, P., Noble, A., Haygarth, P.M., 2023. A global
database of soil plant available phosphorus. Scientific Data 10, 125.
doi:10.1038/s41597-023-02022-4, Table 2. Sample-size-weighted means of
Wolf & Baker (1985), Mallarino (1995) and Buondonno et al. (1992):

    non-calcareous soils:  Olsen = 0.47 * Mehlich3 + 2.4   (n = 337, R2 = 0.80)
    calcareous soils:      Olsen = 0.41 * Mehlich3 + 1.1   (n = 114, R2 = 0.77)

Validated by the authors on 97 USDA-NCSS samples (slope 0.93, NSE > 0.7).

A QUEFTS-specific transfer function exists for tropical soils (Breure,
Van Eynde, Kempen, Comans, Hoffland, 2022. Transfer functions for
phosphorus and potassium soil tests and implications for the QUEFTS model.
Geoderma 406, 115458), ln(P-Olsen) = 0.769 ln(P-M3) + 0.620 ln(Al-M3), but
it requires Mehlich-3 extractable Al, which most soil maps do not provide.
It is therefore not implemented here.
"""

from __future__ import annotations

import numpy as np

#: (slope, intercept) of Olsen P = slope * Mehlich-3 P + intercept, mg/kg (McDowell et al. 2023, Table 2).
MEHLICH3_TO_OLSEN = {
    "non_calcareous": (0.47, 2.4),
    "calcareous": (0.41, 1.1),
}


def mehlich3_to_olsen(p_m3, calcareous=False) -> np.ndarray:
    """Convert Mehlich-3 extractable P to Olsen P (both mg/kg).

    Parameters
    ----------
    p_m3 : Mehlich-3 P, mg/kg (scalar or array; NaN propagates).
    calcareous : bool or boolean array selecting the calcareous equation per
        element. Default False (acid to neutral soils, e.g. most of Honduras).

    Returns
    -------
    Olsen P estimate, mg/kg, never below 0. Note that the intercept gives
    2.4 (or 1.1) mg/kg for a Mehlich-3 value of 0.
    """
    p = np.asarray(p_m3, dtype=float)
    s_nc, i_nc = MEHLICH3_TO_OLSEN["non_calcareous"]
    s_c, i_c = MEHLICH3_TO_OLSEN["calcareous"]
    out = np.where(np.asarray(calcareous, dtype=bool), s_c * p + i_c, s_nc * p + i_nc)
    with np.errstate(invalid="ignore"):
        out = np.where(out < 0, 0.0, out)
    return out


#: Registry of supported P soil-test methods and the converter to Olsen P.
P_METHODS = {
    "olsen": None,
    "mehlich3": mehlich3_to_olsen,
}


def p_to_olsen(p, method: str = "olsen", **kwargs) -> np.ndarray:
    """Convert soil-test P to Olsen P by method name (``"olsen"`` returns the input)."""
    key = method.lower().replace("-", "").replace(" ", "")
    if key not in P_METHODS:
        raise ValueError(f"unknown P method {method!r}; choose from {sorted(P_METHODS)}")
    func = P_METHODS[key]
    p = np.asarray(p, dtype=float)
    return p if func is None else func(p, **kwargs)
