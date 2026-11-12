"""Fertilizer products and application optimization (R-QUEFTS ``ferti_funs.R``).

* :func:`fertilizers` — product table shipped with R-QUEFTS (nutrient % by mass).
* :func:`nutrient_rates` — nutrients supplied by given product amounts.
* :func:`fert_app` — cheapest product mix delivering a nutrient target (linear program).
* :func:`opt_app` — economically optimal application for a QUEFTS model.
"""

from __future__ import annotations

from importlib import resources
from typing import Sequence

import numpy as np
import pandas as pd

from .point import Quefts


def fertilizers() -> pd.DataFrame:
    """Fertilizer product table (``group, name, N, P, K, Ca, Mg, S, Mb, Zn, Bo, Cu`` in % of product)."""
    with resources.files("pyquefts.data").joinpath("fertilizers.csv").open("r", encoding="utf-8") as fh:
        return pd.read_csv(fh)


def nutrient_rates(supply: pd.DataFrame, treatment: Sequence[float]) -> pd.Series:
    """Nutrients (kg) delivered by ``treatment`` kg of each product row of ``supply``.

    ``supply`` must have columns ``N, P, K`` as percentages (R ``nutrientRates``).
    """
    npk = supply[["N", "P", "K"]].to_numpy(dtype=float) * np.asarray(treatment, dtype=float)[:, None] / 100.0
    return pd.Series(npk.sum(axis=0), index=["N", "P", "K"])


def fert_app(nutrients: pd.DataFrame, ferts: pd.DataFrame, price: Sequence[float], exact: bool = True,
             ret_cost: bool = False) -> pd.DataFrame | pd.Series:
    """Least-cost product amounts for target nutrient rates (R ``fertApp``).

    Parameters
    ----------
    nutrients : DataFrame with a row per target and nutrient columns (subset of
        the product columns, e.g. ``N, P, K``), in kg.
    ferts : DataFrame with ``name`` and nutrient percentage columns.
    price : price per kg of each product (len == nrows of ``ferts``).
    exact : if True the nutrient targets must be met exactly (equality
        constraints); if False at least the targets must be supplied.
    ret_cost : return the total cost per target instead of the amounts.

    Returns
    -------
    DataFrame ``name`` + one column per target row with product amounts (kg),
    NaN where no feasible solution exists; or a Series of costs.
    """
    from scipy.optimize import linprog

    price = np.asarray(price, dtype=float)
    if len(price) != len(ferts):
        raise ValueError("price must have one entry per fertilizer product")
    nut_cols = [c for c in nutrients.columns if c in ferts.columns and c != "name"]
    A = ferts[nut_cols].to_numpy(dtype=float).T  # (n_nutrients, n_products), in %
    result = np.full((len(ferts), len(nutrients)), np.nan)
    for i, (_, row) in enumerate(nutrients.iterrows()):
        b = row[nut_cols].to_numpy(dtype=float) * 100.0
        if exact:
            sol = linprog(price, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
        else:
            sol = linprog(price, A_ub=-A, b_ub=-b, bounds=(0, None), method="highs")
        if sol.success:
            result[:, i] = np.where(np.abs(sol.x) < 1e-9, 0.0, sol.x)
    if ret_cost:
        return pd.Series((result * price[:, None]).sum(axis=0), index=nutrients.index)
    out = pd.DataFrame(result, columns=list(nutrients.index))
    out.insert(0, "name", ferts["name"].to_numpy())
    return out


def opt_app(qm: Quefts, ferts: pd.DataFrame, dm_crop_value: float, min_use: float = 0.0,
            max_inv: float = np.inf, start: float = 50.0) -> dict:
    """Economically optimal fertilizer application (R ``optApp``).

    Maximizes net revenue ``(yield - yield0) * dm_crop_value - cost`` over
    product amounts, where yield is the storage-organ yield of the QUEFTS model
    ``qm``. Uses L-BFGS-B with unit finite-difference steps, like R ``optim``.

    Parameters
    ----------
    qm : Quefts model (its fertilizer settings are replaced).
    ferts : DataFrame with ``name, N, P, K`` (as percentages) and ``price_kg``.
    dm_crop_value : value of the product per kg dry matter.
    min_use : amounts below this are set to zero.
    max_inv : maximum total investment.

    Returns
    -------
    dict with ``N, P, K`` (kg/ha nutrients), ``prod0`` (unfertilized yield),
    ``prod`` (yield), ``cost``, ``NR`` (net revenue) and ``fertilizer``
    (DataFrame of selected products and amounts).
    """
    from scipy.optimize import minimize

    need = {"name", "N", "P", "K", "price_kg"}
    if not need.issubset(ferts.columns):
        raise ValueError(f"ferts must have columns {sorted(need)}")
    frt = ferts[["name", "N", "P", "K", "price_kg"]].dropna().reset_index(drop=True)
    if frt.empty:
        raise ValueError("no complete fertilizer records")
    if (frt[["N", "P", "K"]].to_numpy() < 1).all():
        raise ValueError("fertilizer NPK content should be expressed as a percentage")
    content = frt[["N", "P", "K"]].to_numpy(dtype=float) / 100.0
    price = frt["price_kg"].to_numpy(dtype=float)

    def yield_for(p):
        npk = content.T @ p
        qm.fert = qm.fert.replace(N=npk[0], P=npk[1], K=npk[2])
        return qm.run(keep_intermediate=False).scalar("store_lim")

    yield0 = round(yield_for(np.zeros(len(frt))))

    def neg_revenue(p):
        p = np.round(p)
        p[p < min_use] = 0
        cost = float(np.sum(p * price))
        if cost > max_inv:
            return 1.0
        return -((yield_for(p) - yield0) * dm_crop_value - cost)

    n = len(frt)
    res = minimize(neg_revenue, x0=np.full(n, float(start)), method="L-BFGS-B",
                   bounds=[(0, None)] * n, options={"eps": 1.0})
    opt = np.round(res.x)
    opt[opt < min_use] = 0
    npk = content.T @ opt
    qm.fert = qm.fert.replace(N=npk[0], P=npk[1], K=npk[2])
    yld = round(qm.run(keep_intermediate=False).scalar("store_lim"))
    cost = float(np.sum(opt * price))
    nr = (yld - yield0) * dm_crop_value - cost
    sel = pd.DataFrame({"type": frt["name"], "amount": opt, "price_kg": price})
    sel = sel[sel["amount"] > 0].reset_index(drop=True)
    return {"N": npk[0], "P": npk[1], "K": npk[2], "prod0": yield0, "prod": yld,
            "cost": round(cost), "NR": round(nr), "fertilizer": sel}
