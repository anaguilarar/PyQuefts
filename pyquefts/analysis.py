"""Extensions beyond R-QUEFTS: nutrient-yield relationships and response analysis.

Everything here is derived by *repeatedly running the unchanged model*; no new
equations are introduced. The functions provide the outputs requested in the
PyQUEFTS brief that R-QUEFTS leaves to the user (response curves, response
surfaces, uptake-yield envelopes, balanced uptake requirements).
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from .core.model import run_quefts
from .core.parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from .crops import get_crop


def _prep(crop, soil, fert, biom):
    if isinstance(crop, str):
        crop = get_crop(crop)
    return (_coerce(soil, SoilParameters, SoilParameters), _coerce(crop, CropParameters, CropParameters),
            _coerce(fert, Fertilizer, Fertilizer), _coerce(biom, Biomass, Biomass))


def response_curve(nutrient: str, rates: Sequence[float], crop, soil=None, fert=None, biom=None,
                   outputs: Sequence[str] = ("store_lim", "N_uptake", "P_uptake", "K_uptake")) -> pd.DataFrame:
    """Yield (and other outputs) versus the application rate of one nutrient.

    The other two fertilizer rates are held at the values in ``fert``.
    Reproduces the "Explore" example of the reagro.org tutorial.
    """
    soil, crop, fert, biom = _prep(crop, soil, fert, biom)
    nutrient = nutrient.upper()
    if nutrient not in ("N", "P", "K"):
        raise ValueError("nutrient must be 'N', 'P' or 'K'")
    rates = np.asarray(rates, dtype=float)
    fert = fert.replace(**{nutrient: rates})
    res = run_quefts(soil, crop, fert, biom)
    df = pd.DataFrame({nutrient: rates})
    for k in outputs:
        df[k] = res[k]
    return df


def response_surface(n_rates: Sequence[float], p_rates: Sequence[float], k_rates: Sequence[float], crop,
                     soil=None, biom=None, output: str = "store_lim") -> "xarray.DataArray":
    """Full factorial N x P x K response as an ``xarray.DataArray`` (dims N, P, K)."""
    import xarray as xr

    soil, crop, _, biom = _prep(crop, soil, None, biom)
    N, P, K = np.meshgrid(np.asarray(n_rates, float), np.asarray(p_rates, float), np.asarray(k_rates, float),
                          indexing="ij")
    res = run_quefts(soil, crop, Fertilizer(N=N, P=P, K=K), biom, keep_intermediate=False)
    return xr.DataArray(res[output], dims=("N", "P", "K"),
                        coords={"N": np.asarray(n_rates, float), "P": np.asarray(p_rates, float),
                                "K": np.asarray(k_rates, float)}, name=output,
                        attrs={"units": "kg/ha", "crop": crop.name or "default"})


def uptake_yield_envelope(crop, biom=None, uptake_max: Mapping[str, float] | None = None, n: int = 101) -> pd.DataFrame:
    """Maximum-accumulation and maximum-dilution yield lines per nutrient (Janssen 1990 Fig. 1).

    For the given crop and attainable biomass, returns ``Y_acc`` and ``Y_dil``
    (kg/ha of the modelled organ) as linear functions of uptake, using the
    crop-specific ratios ``YratA``/``YratD`` and zero-uptakes that the model
    derives internally. The model itself is run once to obtain those ratios.
    """
    soil, crop, fert, biom = _prep(crop, None, None, biom)
    res = run_quefts(soil, crop, fert, biom)
    frames = []
    for nut in ("N", "P", "K"):
        zero = float(res[f"{nut}zero"])
        ra, rd = float(res[f"{nut}YratA"]), float(res[f"{nut}YratD"])
        umax = (uptake_max or {}).get(nut, float(res[f"{nut}_req_acc"]) * 1.2)
        u = np.linspace(0, umax, n)
        frames.append(pd.DataFrame({"nutrient": nut, "uptake": u,
                                    "Y_acc": ra * np.maximum(0, u - zero),
                                    "Y_dil": rd * np.maximum(0, u - zero)}))
    return pd.concat(frames, ignore_index=True)


def balanced_uptake(crop, target_yields: Sequence[float], ymax: float, leaf_ratio: float, stem_ratio: float,
                    season_length: float = 120.0, soil=None, tol: float = 1e-3) -> pd.DataFrame:
    """Balanced N, P, K uptake required for target yields at a given potential yield.

    Follows Witt et al. (1999) and Yang et al. (2017): QUEFTS is run with the
    attainable (potential) yield ``ymax`` and nutrient supplies that are
    *balanced*, i.e. all three supplies are the same multiple of the mid-range
    requirement. For each target yield the supply multiple is found by
    bisection such that the model yield equals the target; the corresponding
    uptakes are the balanced requirements. The uptake-yield curve is
    linear-parabolic-plateau: linear up to roughly 60-70 % of ``ymax``.

    Returns one row per target with ``N_uptake, P_uptake, K_uptake`` (kg/ha)
    and internal efficiencies ``IE_x`` (kg storage organ per kg uptake).

    Notes
    -----
    * No equation is added: the unchanged model is inverted numerically.
    * For N-fixing crops (``Nfix > 0``) the N uptake is the soil-derived part
      only, as in the model. Pass a crop with ``Nfix=0`` for total N.
    """
    soil, crop, _, _ = _prep(crop, soil, None, None)
    ymax = float(ymax)
    biom = Biomass.from_yield(ymax, leaf_ratio, stem_ratio, SeasonLength=season_length)
    res0 = run_quefts(soil, crop, Fertilizer(0, 0, 0), biom)
    mid = np.array([(res0[f"{n}_req_acc"] + res0[f"{n}_req_dil"]) / 2 for n in "NPK"], dtype=float).ravel()
    fert0 = Fertilizer(0, 0, 0)

    def run_scale(s):
        sl = soil.replace(N_base_supply=mid[0] * s, P_base_supply=mid[1] * s, K_base_supply=mid[2] * s)
        return run_quefts(sl, crop, fert0, biom, keep_intermediate=False)

    rows = []
    for ya in np.asarray(target_yields, dtype=float):
        if ya > ymax:
            rows.append({"target_yield": ya, "yield": np.nan, "N_uptake": np.nan, "P_uptake": np.nan,
                         "K_uptake": np.nan, "IE_N": np.nan, "IE_P": np.nan, "IE_K": np.nan})
            continue
        lo, hi = 0.0, 1.0
        while float(run_scale(hi)["store_lim"]) < ya and hi < 1e4:
            hi *= 2
        for _ in range(80):
            s = (lo + hi) / 2
            if float(run_scale(s)["store_lim"]) >= ya:
                hi = s
            else:
                lo = s
            if hi - lo < tol * max(hi, 1e-9):
                break
        r = run_scale(hi)
        un, up, uk = (float(r[f"{n}_uptake"]) for n in "NPK")
        rows.append({"target_yield": ya, "yield": float(r["store_lim"]), "N_uptake": un, "P_uptake": up,
                     "K_uptake": uk, "IE_N": ya / un if un else np.nan, "IE_P": ya / up if up else np.nan,
                     "IE_K": ya / uk if uk else np.nan})
    return pd.DataFrame(rows)


def limiting_nutrient(result) -> np.ndarray:
    """Name of the nutrient whose dilution yield is lowest (most limiting) per element.

    Uses the yield ranges ``YND, YPD, YKD`` returned with ``keep_intermediate=True``.
    """
    stack = np.stack([np.asarray(result["YND"]), np.asarray(result["YPD"]), np.asarray(result["YKD"])], axis=0)
    idx = np.argmin(stack, axis=0)
    names = np.array(["N", "P", "K"])
    out = names[idx].astype(object)
    out = np.where(np.isnan(stack).any(axis=0), None, out)
    return out
