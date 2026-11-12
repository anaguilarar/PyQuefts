"""Steps 2-4 of QUEFTS: nutrient uptake, yield ranges and nutrient-limited yield.

This is a vectorized (NumPy) transcription of ``QueftsModel::run()`` in
``src/quefts.cpp`` of R-QUEFTS (Hijmans, 2016-2025), which itself implements
Janssen et al. (1990) Tables IV and V with crop-specific yield/uptake ratios
derived from the nutrient concentration limits of vegetative and storage
organs.

Every function accepts scalars or broadcastable arrays. Branches of the
original ``if/else`` logic are evaluated for all elements and selected with
``numpy.where``; divisions that are undefined in an unselected branch are
silenced with ``numpy.errstate``.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from .result import QueftsResult

__all__ = ["uptake", "pair_yield", "requirements", "season_factor", "run_quefts"]


def season_factor(uptake_adjust: np.ndarray, season_length) -> np.ndarray:
    """Multiplier of the base supply for a season of ``season_length`` days.

    Linear interpolation of the (days, factor) table, clamped at both ends.

    Note: the C++ ``approx()`` of R-QUEFTS drops the last table row when
    interpolating. For the default table the last two factors are equal so the
    results are identical; for custom tables PyQUEFTS uses the full table.
    """
    ua = np.asarray(uptake_adjust, dtype=float)
    return np.interp(np.asarray(season_length, dtype=float), ua[:, 0], ua[:, 1])


def uptake(Si, iZero, iYratA, iYratD, Sj, jZero, jYratA, jYratD) -> np.ndarray:
    """Actual uptake of nutrient i given the supply of nutrient j (Janssen 1990, Table IV).

    Parameters
    ----------
    Si, Sj : potential supply of nutrients i and j (kg/ha)
    iZero, jZero : uptake at zero yield (the ``r`` constants of Janssen 1990)
    iYratA, iYratD : yield/uptake ratio of nutrient i at maximum accumulation
        and maximum dilution (``a`` and ``d`` of Janssen 1990); same for j.

    Returns
    -------
    U_i(j): uptake of i limited by the supply of j.
    """
    Si, iZero, iYratA, iYratD, Sj, jZero, jYratA, jYratD = np.broadcast_arrays(
        *(np.asarray(v, dtype=float) for v in (Si, iZero, iYratA, iYratD, Sj, jZero, jYratA, jYratD)))
    dSj = Sj - jZero
    with np.errstate(all="ignore"):
        limA = iZero + jYratA * dSj / iYratD                         # situation A boundary
        limC = iZero + dSj * (2.0 * jYratD / iYratA - jYratA / iYratD)  # situation C boundary
        upB = Si - (0.25 * (Si - iZero - jYratA * dSj / iYratD) ** 2) / (
            (jYratD / iYratA - jYratA / iYratD) * dSj)
        upC = iZero + jYratD * dSj / iYratA
        condA = Si <= limA
        condB = ~condA & (Si < limC)
        up = np.where(condA, Si, np.where(condB, upB, upC))
    return up


def pair_yield(iYA, jYA, iUPT, izero, iYratD, iYratA, iYD, jYD, zYD, YldO) -> np.ndarray:
    """Yield for the nutrient pair (i, j) capped by nutrient z and attainable yield.

    Transcription of ``yield()`` in ``quefts.cpp`` (Janssen 1990, Table V with
    boundary rules). The parabola in the uptake of i runs from ``jYA`` (where
    the dilution yield of i equals the accumulation yield of j) to the upper
    limit ``L = min(iYD, jYD, zYD, YldO)``.

    Parameters
    ----------
    iYA, iYD : yield at max accumulation / dilution of nutrient i
    jYA, jYD : same for nutrient j
    zYD : yield at max dilution of the third nutrient
    iUPT, izero : uptake of i and its zero-yield uptake
    iYratA, iYratD : yield/uptake ratios of i
    YldO : attainable yield of the modelled organ
    """
    iYA, jYA, iUPT, izero, iYratD, iYratA, iYD, jYD, zYD, YldO = np.broadcast_arrays(
        *(np.asarray(v, dtype=float) for v in (iYA, jYA, iUPT, izero, iYratD, iYratA, iYD, jYD, zYD, YldO)))

    LYD = np.minimum(np.minimum(iYD, jYD), np.minimum(zYD, YldO))
    m1 = iYD <= jYA                       # i so diluted that yield is the lower limit
    m2 = ~m1 & (iYA >= LYD)               # i accumulated beyond the upper limit
    m3 = ~m1 & ~m2 & (iYD <= LYD)         # i is itself the most limiting: drop it from L
    LYD2 = np.minimum(jYD, np.minimum(zYD, YldO))
    L = np.where(m3, LYD2, LYD)

    with np.errstate(all="ignore"):
        x = iUPT - izero - jYA / iYratD
        den = L / iYratA - jYA / iYratD
        parab = jYA + 2.0 * (L - jYA) * x / den - (L - jYA) * x ** 2 / den ** 2
    ey = np.where(jYA > L, L, parab)
    ey = np.where(m1 | m2, LYD, ey)
    return ey


def requirements(minVeg, maxVeg, minStore, maxStore, supply, recovery, leaf_att, stem_att, store_att, fix):
    """Nutrient requirements at balanced, accumulated and diluted concentrations.

    Returns ``(reqAcc, reqDil, reqFert)``: uptake needed for the attainable
    biomass at maximum accumulation, at maximum dilution, and the fertilizer
    needed at balanced (mid-range) concentrations given the supply and
    recovery fraction. ``fix`` is the N-fixation fraction (0 for P and K);
    following R-QUEFTS the fertilizer requirement is zero whenever ``fix != 0``.
    """
    (minVeg, maxVeg, minStore, maxStore, supply, recovery, leaf_att, stem_att, store_att, fix) = (
        np.asarray(v, dtype=float) for v in
        (minVeg, maxVeg, minStore, maxStore, supply, recovery, leaf_att, stem_att, store_att, fix))
    veg = leaf_att + stem_att
    limVeg = (minVeg + maxVeg) / 2.0
    limStore = (maxStore + minStore) / 2.0
    req = np.where(fix != 0, 0.0, limVeg * veg + limStore * store_att)
    with np.errstate(all="ignore"):
        reqFert = np.maximum(0.0, (req - supply) / recovery)
    reqAcc = (maxVeg * veg + maxStore * store_att) * (1.0 - fix)
    reqDil = (minVeg * veg + minStore * store_att) * (1.0 - fix)
    return reqAcc, reqDil, reqFert


def run_quefts(soil: SoilParameters | Mapping | None = None, crop: CropParameters | Mapping | None = None,
               fert: Fertilizer | Mapping | None = None, biom: Biomass | Mapping | None = None,
               keep_intermediate: bool = True) -> QueftsResult:
    """Run the QUEFTS model (steps 2-4) for scalar or array inputs.

    Parameters
    ----------
    soil, crop, fert, biom : parameter objects or mappings with the R-QUEFTS
        field names. ``None`` uses the R-QUEFTS defaults.
    keep_intermediate : also return pairwise uptakes, yield ranges and the
        six pairwise yields (needed for nutrient-yield relationship outputs).

    Returns
    -------
    QueftsResult with the 12 R-QUEFTS outputs (``leaf_lim, stem_lim,
    store_lim, N_supply, P_supply, K_supply, N_uptake, P_uptake, K_uptake,
    N_gap, P_gap, K_gap``) and, optionally, the intermediates.

    Elements with NaN in any supply or biomass input yield NaN outputs.
    """
    soil = _coerce(soil, SoilParameters, SoilParameters)
    crop = _coerce(crop, CropParameters, CropParameters)
    fert = _coerce(fert, Fertilizer, Fertilizer)
    biom = _coerce(biom, Biomass, Biomass)

    f = lambda v: np.asarray(v, dtype=float)  # noqa: E731
    leaf_att, stem_att, store_att = f(biom.leaf_att), f(biom.stem_att), f(biom.store_att)
    Yzero, Nfix = f(crop.Yzero), f(crop.Nfix)

    # ---- which organ is modelled
    organs3 = store_att > 200.0
    organs1 = ~organs3 & (stem_att <= 200.0)
    organs2 = ~organs3 & ~organs1
    organ = np.where(organs3, 3, np.where(organs1, 1, 2))
    Yo = np.where(organs3, store_att, np.where(organs1, np.maximum(5.0, leaf_att), stem_att))
    Nzero = np.where(organs3, Yzero * (1.0 - Nfix) * (f(crop.NmaxVeg) + 2.0 * f(crop.NminVeg)) / 3.0, 0.0)
    Pzero = np.where(organs3, Yzero * (f(crop.PmaxVeg) + 2.0 * f(crop.PminVeg)) / 3.0, 0.0)
    Kzero = np.where(organs3, Yzero * (f(crop.KmaxVeg) + 2.0 * f(crop.KminVeg)) / 3.0, 0.0)

    # ---- supply during the season (soil x season factor + fertilizer x recovery)
    relgrow = season_factor(soil.UptakeAdjust, biom.SeasonLength)
    N_supply = f(soil.N_base_supply) * relgrow + f(fert.N) * f(soil.N_recovery)
    P_supply = f(soil.P_base_supply) * relgrow + f(fert.P) * f(soil.P_recovery)
    K_supply = f(soil.K_base_supply) * relgrow + f(fert.K) * f(soil.K_recovery)

    # ---- requirements and fertilizer gaps
    NreqA, NreqD, N_gap = requirements(crop.NminVeg, crop.NmaxVeg, crop.NminStore, crop.NmaxStore, N_supply,
                                       soil.N_recovery, leaf_att, stem_att, store_att, Nfix)
    PreqA, PreqD, P_gap = requirements(crop.PminVeg, crop.PmaxVeg, crop.PminStore, crop.PmaxStore, P_supply,
                                       soil.P_recovery, leaf_att, stem_att, store_att, 0.0)
    KreqA, KreqD, K_gap = requirements(crop.KminVeg, crop.KmaxVeg, crop.KminStore, crop.KmaxStore, K_supply,
                                       soil.K_recovery, leaf_att, stem_att, store_att, 0.0)

    # ---- yield / uptake ratios at accumulation (A) and dilution (D)
    NYratA = Yo / np.maximum(0.1, NreqA - Nzero)
    NYratD = Yo / np.maximum(0.1, NreqD - Nzero)
    PYratA = Yo / np.maximum(0.1, PreqA - Pzero)
    PYratD = Yo / np.maximum(0.1, PreqD - Pzero)
    KYratA = Yo / np.maximum(0.1, KreqA - Kzero)
    KYratD = Yo / np.maximum(0.1, KreqD - Kzero)

    # ---- step 2: uptake
    UNP = uptake(N_supply, Nzero, NYratA, NYratD, P_supply, Pzero, PYratA, PYratD)
    UNK = uptake(N_supply, Nzero, NYratA, NYratD, K_supply, Kzero, KYratA, KYratD)
    UN = np.maximum(0.0, np.minimum(UNP, UNK))
    UPN = uptake(P_supply, Pzero, PYratA, PYratD, N_supply, Nzero, NYratA, NYratD)
    UPK = uptake(P_supply, Pzero, PYratA, PYratD, K_supply, Kzero, KYratA, KYratD)
    UP = np.maximum(0.0, np.minimum(UPN, UPK))
    UKN = uptake(K_supply, Kzero, KYratA, KYratD, N_supply, Nzero, NYratA, NYratD)
    UKP = uptake(K_supply, Kzero, KYratA, KYratD, P_supply, Pzero, PYratA, PYratD)
    UK = np.maximum(0.0, np.minimum(UKN, UKP))

    # ---- step 3: yield ranges
    YNA = NYratA * np.maximum(0.0, UN - Nzero)
    YND = NYratD * np.maximum(0.0, UN - Nzero)
    YPA = PYratA * np.maximum(0.0, UP - Pzero)
    YPD = PYratD * np.maximum(0.0, UP - Pzero)
    YKA = KYratA * np.maximum(0.0, UK - Kzero)
    YKD = KYratD * np.maximum(0.0, UK - Kzero)

    # ---- step 4: pairwise yields and final yield
    EYNP = pair_yield(YNA, YPA, UN, Nzero, NYratD, NYratA, YND, YPD, YKD, Yo)
    EYPN = pair_yield(YPA, YNA, UP, Pzero, PYratD, PYratA, YPD, YND, YKD, Yo)
    EYNK = pair_yield(YNA, YKA, UN, Nzero, NYratD, NYratA, YND, YKD, YPD, Yo)
    EYKN = pair_yield(YKA, YNA, UK, Kzero, KYratD, KYratA, YKD, YND, YPD, Yo)
    EYPK = pair_yield(YPA, YKA, UP, Pzero, PYratD, PYratA, YPD, YKD, YND, Yo)
    EYKP = pair_yield(YKA, YPA, UK, Kzero, KYratD, KYratA, YKD, YPD, YND, Yo)

    limYield = np.minimum(np.minimum(YND, YPD), np.minimum(YKD, Yo))
    yld = np.minimum(limYield, (EYNP + EYPN + EYNK + EYKN + EYPK + EYKP) / 6.0)

    # ---- partition over organs
    with np.errstate(all="ignore"):
        veg = leaf_att + stem_att
        leaf3 = leaf_att / veg * ((veg - Yzero) * yld / store_att + Yzero)
        # NOTE: R-QUEFTS uses (NminVeg + 2*NminVeg) here; (NmaxVeg + 2*NminVeg) was
        # probably intended (it would match Nzero). Reproduced verbatim for fidelity.
        leaf3_zero = leaf_att / veg * np.minimum.reduce([
            np.broadcast_to(Yzero, yld.shape).astype(float),
            UN / (0.333 * (f(crop.NminVeg) + 2.0 * f(crop.NminVeg)) * (1.0 - Nfix)),
            UP / (0.333 * (f(crop.PminVeg) + 2.0 * f(crop.PminVeg))),
            UK / (0.333 * (f(crop.KminVeg) + 2.0 * f(crop.KminVeg)))])
        leaf3 = np.where(yld <= 0, leaf3_zero, leaf3)
        stem3 = (stem_att / leaf_att) * leaf3

        leaf1 = yld
        stem1 = (stem_att / leaf_att) * leaf1
        store1 = (store_att / leaf_att) * leaf1

        stem2 = yld
        leaf2 = (leaf_att / stem_att) * stem2
        store2 = (store_att / stem_att) * stem2

    leaf_lim = np.where(organs3, leaf3, np.where(organs1, leaf1, leaf2))
    stem_lim = np.where(organs3, stem3, np.where(organs1, stem1, stem2))
    store_lim = np.where(organs3, yld, np.where(organs1, store1, store2))

    out = {
        "leaf_lim": leaf_lim, "stem_lim": stem_lim, "store_lim": store_lim,
        "N_supply": N_supply, "P_supply": P_supply, "K_supply": K_supply,
        "N_uptake": UN, "P_uptake": UP, "K_uptake": UK,
        "N_gap": N_gap, "P_gap": P_gap, "K_gap": K_gap,
    }
    if keep_intermediate:
        out.update({
            "season_factor": relgrow, "organ": organ, "Y_att": Yo,
            "Nzero": Nzero, "Pzero": Pzero, "Kzero": Kzero,
            "N_req_acc": NreqA, "N_req_dil": NreqD, "P_req_acc": PreqA, "P_req_dil": PreqD,
            "K_req_acc": KreqA, "K_req_dil": KreqD,
            "NYratA": NYratA, "NYratD": NYratD, "PYratA": PYratA, "PYratD": PYratD,
            "KYratA": KYratA, "KYratD": KYratD,
            "N_uptake_P": UNP, "N_uptake_K": UNK, "P_uptake_N": UPN, "P_uptake_K": UPK,
            "K_uptake_N": UKN, "K_uptake_P": UKP,
            "YNA": YNA, "YND": YND, "YPA": YPA, "YPD": YPD, "YKA": YKA, "YKD": YKD,
            "Y_NP": EYNP, "Y_PN": EYPN, "Y_NK": EYNK, "Y_KN": EYKN, "Y_PK": EYPK, "Y_KP": EYKP,
            "yield_lim": yld,
        })

    # ---- propagate missing inputs
    invalid = np.isnan(N_supply) | np.isnan(P_supply) | np.isnan(K_supply) | np.isnan(Yo)
    invalid = np.broadcast_to(invalid, yld.shape)
    if invalid.any():
        for k, v in out.items():
            v = np.broadcast_to(np.asarray(v, dtype=float), yld.shape)
            out[k] = np.where(invalid, np.nan, v)
    else:
        for k, v in out.items():
            out[k] = np.broadcast_to(np.asarray(v), yld.shape).copy()

    return QueftsResult(out)
