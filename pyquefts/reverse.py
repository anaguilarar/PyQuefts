"""Reverse QUEFTS: estimate apparent soil nutrient supply from omission-trial yields (R ``revSupply``)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .core.model import run_quefts
from .core.parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from .crops import get_crop


def rev_supply(obs: pd.DataFrame, crop, soil=None, Ya: float = None, leaf_ratio: float = None,
               stem_ratio: float = None, season_length: float = 120.0, **optim_kwargs) -> np.ndarray:
    """Fit N, P, K base supplies that minimize the squared yield error of a trial.

    Parameters
    ----------
    obs : DataFrame with columns ``N, P, K`` (fertilizer, kg/ha) and ``Y``
        (observed storage-organ yield, kg/ha).
    crop : CropParameters, mapping or crop name.
    soil : SoilParameters giving recoveries, UptakeAdjust and the initial
        guess of the base supplies.
    Ya : attainable yield (kg/ha).
    leaf_ratio, stem_ratio : leaf and stem biomass relative to ``Ya``.
    season_length : days.
    optim_kwargs : passed to ``scipy.optimize.minimize`` (method L-BFGS-B).

    Returns
    -------
    ndarray ``[N_base_supply, P_base_supply, K_base_supply]`` (kg/ha).
    """
    from scipy.optimize import minimize

    if isinstance(crop, str):
        crop = get_crop(crop)
    crop = _coerce(crop, CropParameters, CropParameters)
    soil = _coerce(soil, SoilParameters, SoilParameters)
    obs = obs[["N", "P", "K", "Y"]].astype(float)
    fert = Fertilizer(N=obs["N"].to_numpy(), P=obs["P"].to_numpy(), K=obs["K"].to_numpy())
    biom = Biomass.from_yield(Ya, leaf_ratio, stem_ratio, SeasonLength=season_length)
    y = obs["Y"].to_numpy()

    def sse(s):
        sl = soil.replace(N_base_supply=s[0], P_base_supply=s[1], K_base_supply=s[2])
        yq = run_quefts(sl, crop, fert, biom, keep_intermediate=False)["store_lim"]
        return float(np.sum((yq - y) ** 2))

    x0 = np.array([soil.N_base_supply, soil.P_base_supply, soil.K_base_supply], dtype=float)
    res = minimize(sse, x0, method="L-BFGS-B", bounds=[(0, None)] * 3, **optim_kwargs)
    return res.x
