"""Step 1 of QUEFTS: potential soil nutrient supply from topsoil properties.

Two parameterizations are provided, mirroring R-QUEFTS:

* :func:`nut_supply1` — Janssen et al. (1990), Geoderma 46:299-318, Table II.
* :func:`nut_supply2` — Sattari et al. (2014), Field Crops Research 157:35-46,
  Table 4, as implemented in R-QUEFTS (temperature-dependent N mineralization,
  piecewise pH factors as interpolation tables, total P required).

Units (0-20 cm topsoil): pH in water; SOC in g/kg; exchangeable K in mmol/kg;
P-Olsen and total P in mg/kg; temperature in degrees C. Output: kg/ha of N, P
and K for a 120-day season.

All inputs may be scalars or broadcastable arrays; NaN propagates.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .parameters import ArrayLike, SoilParameters

# --- Sattari (2014) pH factors as used by R-QUEFTS (piecewise linear, clamped) ---
_FN_TABLE = (np.array([4.7, 7.0]), np.array([0.4, 1.0]))
_FP_TABLE = (
    np.array([4.7, 5.0, 5.5, 6.0, 6.7, 7.0, 7.5, 8.0]),
    np.array([0.02, 0.5, 0.875, 1.0, 1.0, 0.9775, 0.84, 0.02]),
)
_FK_TABLE = (np.array([4.5, 5.0, 6.0, 6.8]), np.array([1.0, 0.8842312, 0.71, 0.6]))

#: Suggested total-P / P-Olsen ratios for soils where total P was not measured.
PTOTAL_RATIOS = {
    "sattari2014": 95.0,  # Sattari et al. 2014, Kenyan soils without fertilizer P
    "afsis": 55.0,        # Rquefts documentation (AfSIS data)
}


@dataclass(frozen=True, eq=False)
class NutrientSupply:
    """Potential N, P, K supply of the unfertilized soil (kg/ha, 120-day season)."""

    N_base_supply: ArrayLike
    P_base_supply: ArrayLike
    K_base_supply: ArrayLike

    def to_soil(self, **kwargs) -> SoilParameters:
        """Convert to :class:`SoilParameters`; ``kwargs`` set recoveries/UptakeAdjust."""
        return SoilParameters(N_base_supply=self.N_base_supply, P_base_supply=self.P_base_supply,
                              K_base_supply=self.K_base_supply, **kwargs)

    def to_dict(self) -> dict:
        return {"N_base_supply": self.N_base_supply, "P_base_supply": self.P_base_supply,
                "K_base_supply": self.K_base_supply}

    def to_array(self) -> np.ndarray:
        """Stack as (..., 3) array in N, P, K order."""
        return np.stack(np.broadcast_arrays(self.N_base_supply, self.P_base_supply, self.K_base_supply), axis=-1)

    def __getitem__(self, key):
        return getattr(self, key)


# ---------------------------------------------------------------- Janssen 1990

def ph_factors_janssen(ph: ArrayLike) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """pH correction factors fN, fP, fK of Janssen et al. (1990) eqs. 4-6."""
    ph = np.asarray(ph, dtype=float)
    fN = 0.25 * (ph - 3.0)
    fP = 1.0 - 0.5 * (ph - 6.0) ** 2
    fK = 0.625 * (3.4 - 0.4 * ph)
    return fN, fP, fK


def nut_supply1(ph: ArrayLike, soc: ArrayLike, kex: ArrayLike, polsen: ArrayLike,
                ptotal: Optional[ArrayLike] = None) -> NutrientSupply:
    """Potential supply after Janssen et al. (1990), Table II (R ``nutSupply1``).

    Parameters
    ----------
    ph : soil pH (H2O)
    soc : soil organic carbon, g/kg
    kex : exchangeable K, mmol/kg
    polsen : P-Olsen, mg/kg
    ptotal : total P, mg/kg; optional. Where given (not NaN) eq. 8 is used
        instead of eq. 2.

    Notes
    -----
    Negative supplies are clamped to zero, as in R-QUEFTS.
    """
    ph, soc, kex, polsen = (np.asarray(v, dtype=float) for v in (ph, soc, kex, polsen))
    fN, fP, fK = ph_factors_janssen(ph)

    n_sup = fN * 6.8 * soc                                   # eq. 1
    p_sup = fP * 0.35 * soc + 0.5 * polsen                   # eq. 2
    if ptotal is not None:
        ptotal = np.asarray(ptotal, dtype=float)
        p_tot = fP * 0.014 * ptotal + 0.5 * polsen           # eq. 8
        p_sup = np.where(np.isnan(ptotal), p_sup, p_tot)
    k_sup = (fK * 400.0 * kex) / (2.0 + 0.9 * soc)           # eq. 3

    n_sup, p_sup, k_sup = np.broadcast_arrays(n_sup, p_sup, k_sup)
    with np.errstate(invalid="ignore"):
        n_sup = np.where(n_sup < 0, 0.0, n_sup)
        p_sup = np.where(p_sup < 0, 0.0, p_sup)
        k_sup = np.where(k_sup < 0, 0.0, k_sup)
    return NutrientSupply(n_sup, p_sup, k_sup)


# ---------------------------------------------------------------- Sattari 2014

def ph_factors_sattari(ph: ArrayLike) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """pH factors fN, fP, fK as interpolation tables (R-QUEFTS ``.get_fn/.get_fp/.get_fk``).

    These are the piecewise-linear approximations R-QUEFTS uses for the
    piecewise functions of Sattari et al. (2014), Table 4. ``numpy.interp``
    clamps outside the table range, equivalent to ``approxfun(rule=2)``.
    """
    ph = np.asarray(ph, dtype=float)
    fN = np.interp(ph, *_FN_TABLE)
    fP = np.interp(ph, *_FP_TABLE)
    fK = np.interp(ph, *_FK_TABLE)
    return fN, fP, fK


def temperature_factor_n(temp: ArrayLike) -> np.ndarray:
    """Temperature-dependent N mineralization coefficient ``2 * 2^((T-9)/9)`` (Sattari 2014, Table 4).

    Replaces the constant 6.8 of Janssen (1990); equals 6.8 at T ~ 24.9 C.
    """
    temp = np.asarray(temp, dtype=float)
    return 2.0 * 2.0 ** ((temp - 9.0) / 9.0)


def nut_supply2(temp: ArrayLike, ph: ArrayLike, soc: ArrayLike, kex: ArrayLike, polsen: ArrayLike,
                ptotal: Optional[ArrayLike] = None, ptotal_ratio: Optional[float | str] = None) -> NutrientSupply:
    """Potential supply after Sattari et al. (2014) as implemented in R ``nutSupply2``.

    Parameters
    ----------
    temp : mean growing-season temperature, C
    
    ph : soil pH (H2O)
    soc : soil organic carbon, g/kg
    kex : exchangeable K, mmol/kg
    polsen : P-Olsen, mg/kg
    ptotal : total P, mg/kg. Required by the method. If None, it is estimated
        as ``ptotal_ratio * polsen``.
    ptotal_ratio : float or key of :data:`PTOTAL_RATIOS` ("sattari2014" = 95,
        "afsis" = 55). Only used when ``ptotal`` is None. Defaults to "afsis",
        the value recommended in the Rquefts documentation.

    Notes
    -----
    No clamping of negative values is applied (same as R-QUEFTS).
    """
    temp, ph, soc, kex, polsen = (np.asarray(v, dtype=float) for v in (temp, ph, soc, kex, polsen))
    if ptotal is None:
        if ptotal_ratio is None:
            ptotal_ratio = "afsis"
        if isinstance(ptotal_ratio, str):
            ptotal_ratio = PTOTAL_RATIOS[ptotal_ratio]
        ptotal = float(ptotal_ratio) * polsen
    ptotal = np.asarray(ptotal, dtype=float)

    fN, fP, fK = ph_factors_sattari(ph)
    n_sup = fN * temperature_factor_n(temp) * soc
    p_sup = fP * 0.014 * ptotal + 0.5 * polsen
    k_sup = (fK * 400.0 * kex) / (2.0 + 0.9 * soc)
    n_sup, p_sup, k_sup = np.broadcast_arrays(n_sup, p_sup, k_sup)
    return NutrientSupply(np.array(n_sup), np.array(p_sup), np.array(k_sup))


def nut_supply(method: str, **kwargs) -> NutrientSupply:
    """Dispatch to :func:`nut_supply1` (``"janssen1990"``) or :func:`nut_supply2` (``"sattari2014"``)."""
    method = method.lower()
    if method in ("1", "janssen", "janssen1990", "nutsupply1"):
        return nut_supply1(**kwargs)
    if method in ("2", "sattari", "sattari2014", "nutsupply2"):
        return nut_supply2(**kwargs)
    raise ValueError(f"unknown supply method {method!r}")
