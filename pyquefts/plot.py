"""Plotting helpers (extension; requires matplotlib, ``pip install pyquefts[plot]``).

Nothing here changes the model: every function plots the outputs of
:func:`pyquefts.core.model.run_quefts` (or of the :mod:`pyquefts.analysis`
helpers, which only call it repeatedly) against the inputs that were varied,
the way the Rquefts tutorial does with base-R ``plot()``.

Typical use::

    q = Quefts(soil=sup.to_soil(), crop="maize", fert=dict(N=np.arange(0, 160, 10), P=10, K=60),
               biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
    q.run().plot()                          # yield vs the input that varies (N here)
    q.run().plot(y=["store_lim", "N_uptake", "N_gap"])
    plot_response_curves(q, rates=np.arange(0, 201, 10))      # N, P and K panels
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .core.model import run_quefts
from .core.parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from .core.result import R_OUTPUTS, QueftsResult
from .crops import get_crop

__all__ = ["plot_result", "plot_response_curves", "plot_curves", "plot_spatial_response", "plot_surface",
           "plot_envelope", "LABELS"]

#: Axis labels for inputs and outputs.
LABELS = {
    "N": "N applied (kg/ha)", "P": "P applied (kg/ha)", "K": "K applied (kg/ha)",
    "store_att": "attainable storage-organ yield (kg/ha DM)", "leaf_att": "attainable leaf biomass (kg/ha DM)",
    "stem_att": "attainable stem biomass (kg/ha DM)", "SeasonLength": "season length (days)",
    "N_base_supply": "soil N supply (kg/ha)", "P_base_supply": "soil P supply (kg/ha)",
    "K_base_supply": "soil K supply (kg/ha)",
    "N_recovery": "N fertilizer recovery (-)", "P_recovery": "P fertilizer recovery (-)",
    "K_recovery": "K fertilizer recovery (-)",
    "store_lim": "storage-organ yield (kg/ha DM)", "leaf_lim": "leaf biomass (kg/ha DM)",
    "stem_lim": "stem biomass (kg/ha DM)",
    "N_supply": "N supply (kg/ha)", "P_supply": "P supply (kg/ha)", "K_supply": "K supply (kg/ha)",
    "N_uptake": "N uptake (kg/ha)", "P_uptake": "P uptake (kg/ha)", "K_uptake": "K uptake (kg/ha)",
    "N_gap": "N gap (kg/ha)", "P_gap": "P gap (kg/ha)", "K_gap": "K gap (kg/ha)",
    "uptake": "uptake (kg/ha)",
}

_ORGAN_OUTPUT = {1: "leaf_lim", 2: "stem_lim", 3: "store_lim"}


def _label(key: str) -> str:
    return LABELS.get(key, key)


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise ImportError("plotting requires matplotlib: pip install 'pyquefts[plot]'") from e
    return plt


def _full(result: QueftsResult, values) -> np.ndarray:
    """Broadcast an input or output to the result shape and flatten it."""
    return np.broadcast_to(np.asarray(values, dtype=float), result.shape).ravel()


def _varying_input(result: QueftsResult, x: str | None):
    """Name and values of the x variable: the single input that varies unless ``x`` is given."""
    inputs = getattr(result, "inputs", None) or {}
    if x is not None:
        if x in inputs:
            return x, _full(result, inputs[x])
        if x in result:
            return x, _full(result, result[x])
        raise KeyError(f"{x!r} is neither a recorded input {sorted(inputs)} nor an output")
    cand = [k for k, v in inputs.items() if np.ndim(v) > 0 and np.size(v) > 1]
    if not cand:
        raise ValueError("no input varies in this run (or inputs were not recorded); pass x='<name>'")
    if len(cand) == 1:
        k = cand[0]
    elif set(cand) <= {"leaf_att", "stem_att", "store_att"}:
        k = "store_att"  # biomass built with Biomass.from_yield: all three vary together
    else:
        raise ValueError(f"several inputs vary ({cand}); pass x='<name>' and optionally hue='<name>'")
    return k, _full(result, inputs[k])


def plot_result(result: QueftsResult, x: str | None = None, y: str | Sequence[str] = "store_lim", *,
                hue: str | None = None, ax=None, marker: str = "o", **kwargs):
    """Plot outputs of a vectorized run against the input that was varied.

    Parameters
    ----------
    result : :class:`QueftsResult` from :meth:`pyquefts.Quefts.run` (which records
        the inputs) or any result when ``x`` names an output.
    x : input (``"N"``, ``"P"``, ``"K"``, ``"store_att"``, ``"N_base_supply"``, ...)
        or output to put on the x axis; detected automatically when exactly one
        input varies.
    y : output name or list of names; one panel per output.
    hue : a second varying input; one line per distinct value (e.g. ``x="N", hue="store_att"``).
    ax : matplotlib Axes (or array of Axes, one per ``y``); created when None.

    Returns the Axes (or array of Axes when several ``y``).
    """
    plt = _plt()
    xname, xv = _varying_input(result, x)
    ys = [y] if isinstance(y, str) else list(y)
    if ax is None:
        _, axes = plt.subplots(1, len(ys), figsize=(4.6 * len(ys), 3.8), squeeze=False)
        axes = axes[0]
    else:
        axes = np.atleast_1d(ax).ravel()
        if axes.size != len(ys):
            raise ValueError(f"{len(ys)} outputs requested but {axes.size} axes given")
    groups = [(None, np.ones(xv.size, bool))]
    if hue is not None:
        inputs = getattr(result, "inputs", None) or {}
        hv = _full(result, inputs[hue] if hue in inputs else result[hue])
        groups = [(v, hv == v) for v in np.unique(hv[~np.isnan(hv)])]
    for a, k in zip(axes, ys):
        yv = _full(result, result[k])
        for val, m in groups:
            order = np.argsort(xv[m], kind="stable")
            lbl = None if val is None else f"{hue} = {val:g}"
            a.plot(xv[m][order], yv[m][order], marker=marker, label=lbl, **kwargs)
        a.set_xlabel(_label(xname))
        a.set_ylabel(_label(k))
        a.grid(alpha=0.3)
        if hue is not None:
            a.legend(fontsize="small")
    return axes[0] if len(ys) == 1 else axes


def plot_response_curves(model=None, rates=None, *, nutrients: str | Sequence[str] = "NPK",
                         y: str = "store_lim", crop=None, soil=None, fert=None, biom=None, axes=None,
                         sharey: bool = True, **kwargs):
    """Yield (or another output) against the rate of each nutrient, one panel per nutrient.

    For each nutrient in ``nutrients`` the model is run over ``rates`` with
    the other two fertilizer rates held at their values in ``model.fert``
    (or ``fert``). This reproduces the loop of the Rquefts tutorial
    (``for (i in 1:length(N)) { q["N"] <- N[i]; run(q) }``) in one vectorized
    call per nutrient.

    Parameters
    ----------
    model : :class:`pyquefts.Quefts` instance (soil, crop, fert, biom taken from it);
        alternatively pass ``crop``, ``soil``, ``fert``, ``biom`` directly.
    rates : array of application rates (kg/ha) used for every nutrient, or a
        mapping ``{"N": [...], "P": [...], "K": [...]}``.
    nutrients : subset / order of ``"NPK"``.
    y : output to plot.
    axes : sequence of Axes (one per nutrient); created when None.

    Returns the array of Axes.
    """
    plt = _plt()
    if model is not None:
        soil_p, crop_p, fert_p, biom_p = model.soil, model.crop, model.fert, model.biom
    else:
        if isinstance(crop, str):
            crop = get_crop(crop)
        soil_p = _coerce(soil, SoilParameters, SoilParameters)
        crop_p = _coerce(crop, CropParameters, CropParameters)
        fert_p = _coerce(fert, Fertilizer, Fertilizer)
        biom_p = _coerce(biom, Biomass, Biomass)
    if rates is None:
        rates = np.arange(0, 201, 10)
    nuts = [n.upper() for n in nutrients]
    rate_map = {n: np.asarray(rates[n] if isinstance(rates, Mapping) else rates, dtype=float) for n in nuts}
    if axes is None:
        _, axes = plt.subplots(1, len(nuts), figsize=(4.6 * len(nuts), 3.8), sharey=sharey, squeeze=False)
        axes = axes[0]
    axes = np.atleast_1d(axes).ravel()
    marker = kwargs.pop("marker", "o")
    for a, n in zip(axes, nuts):
        r = rate_map[n]
        base = {}
        for m in "NPK":
            if m == n:
                continue
            v = np.asarray(getattr(fert_p, m), dtype=float)
            if v.size != 1:
                raise ValueError(f"fertilizer {m} must be a scalar to draw the {n} response (it has shape {v.shape}); "
                                 f"set it first, e.g. model.set({m}=<rate>)")
            base[m] = float(v.ravel()[0])
        res = run_quefts(soil_p, crop_p, fert_p.replace(**{n: r}, **base), biom_p, keep_intermediate=y not in R_OUTPUTS)
        a.plot(r, np.asarray(res[y], dtype=float), marker=marker, **kwargs)
        a.set_xlabel(_label(n))
        a.set_ylabel(_label(y))
        a.set_title(", ".join(f"{m} = {v:g}" for m, v in base.items()), fontsize="small")
        a.grid(alpha=0.3)
    name = crop_p.name or "default crop"
    axes[0].figure.suptitle(f"{name}: response to fertilizer (QUEFTS)", fontsize="medium")
    return axes


def plot_curves(df, x: str | None = None, y: str | Sequence[str] = "store_lim", *, hue: str | None = "label",
                ax=None, marker: str = "o", **kwargs):
    """Plot a long response table (one line per ``hue`` group), e.g. from :func:`pyquefts.geo.response_points`.

    ``x`` defaults to ``df.attrs["nutrient"]`` (set by ``response_points``);
    ``hue`` is a column naming the groups (``"label"`` = "lon, lat"; use
    ``hue=None`` for a single line). One panel per entry of ``y``.
    """
    plt = _plt()
    if x is None:
        x = df.attrs.get("nutrient")
        if x is None:
            raise ValueError("pass x='<column>' (the table has no 'nutrient' attribute)")
    ys = [y] if isinstance(y, str) else list(y)
    if ax is None:
        _, axes = plt.subplots(1, len(ys), figsize=(4.6 * len(ys), 3.8), squeeze=False)
        axes = axes[0]
    else:
        axes = np.atleast_1d(ax).ravel()
    groups = [(None, df)] if hue is None or hue not in df else list(df.groupby(hue, sort=False))
    for a, k in zip(axes, ys):
        for lbl, g in groups:
            g = g.sort_values(x)
            a.plot(g[x].to_numpy(), g[k].to_numpy(), marker=marker, label=None if lbl is None else str(lbl), **kwargs)
        a.set_xlabel(_label(x))
        a.set_ylabel(_label(k))
        a.grid(alpha=0.3)
        if hue is not None and len(groups) > 1:
            a.legend(fontsize="small", title=hue)
    return axes[0] if len(ys) == 1 else axes


def plot_spatial_response(ds, y: str = "store_lim", nutrient: str | None = None, *, quantiles=(0.1, 0.9), ax=None,
                          **kwargs):
    """Area-mean response curve with a quantile band from :func:`pyquefts.geo.response_spatial`.

    ``ds`` has a rate dimension named after the nutrient (detected when
    ``nutrient`` is None); the mean and the ``quantiles`` of ``y`` are taken
    over all grid cells for each rate.
    """
    plt = _plt()
    if nutrient is None:
        cand = [d for d in ds[y].dims if d in ("N", "P", "K")]
        if len(cand) != 1:
            raise ValueError("pass nutrient='N'|'P'|'K' (rate dimension not found)")
        nutrient = cand[0]
    da = ds[y]
    space = [d for d in da.dims if d != nutrient]
    stacked = da.stack(cell=space)
    mean = stacked.mean("cell", skipna=True).compute()
    lo = stacked.quantile(quantiles[0], dim="cell", skipna=True).compute()
    hi = stacked.quantile(quantiles[1], dim="cell", skipna=True).compute()
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 3.8))
    r = ds[nutrient].values
    ax.fill_between(r, lo.values, hi.values, alpha=0.2, label=f"{quantiles[0]:g}-{quantiles[1]:g} quantiles")
    ax.plot(r, mean.values, marker=kwargs.pop("marker", "o"), label="area mean", **kwargs)
    ax.set_xlabel(_label(nutrient))
    ax.set_ylabel(_label(y))
    ax.grid(alpha=0.3)
    ax.legend(fontsize="small")
    return ax


def plot_surface(da, x: str = "N", y: str = "P", fixed: float | None = None, *, ax=None, levels=15,
                 cmap: str = "viridis", contour_lines: bool = True, **kwargs):
    """Filled contours of a :func:`pyquefts.analysis.response_surface` result.

    ``da`` has dims ``N, P, K``; the dimension that is neither ``x`` nor ``y``
    is fixed at ``fixed`` (nearest value) or at its first value.
    """
    plt = _plt()
    other = [d for d in da.dims if d not in (x, y)]
    fixed_txt = ""
    if other:
        d = other[0]
        sel = da[d].values[0] if fixed is None else fixed
        da = da.sel({d: sel}, method="nearest")
        fixed_txt = f"{d} = {float(da[d]):g} kg/ha"
    if ax is None:
        _, ax = plt.subplots(figsize=(5.2, 4.2))
    z = da.transpose(y, x).values
    cf = ax.contourf(da[x].values, da[y].values, z, levels=levels, cmap=cmap, **kwargs)
    if contour_lines:
        cs = ax.contour(da[x].values, da[y].values, z, levels=cf.levels, colors="k", linewidths=0.5, alpha=0.6)
        ax.clabel(cs, fontsize="x-small", fmt="%g")
    ax.figure.colorbar(cf, ax=ax, label=_label(da.name or "value"))
    ax.set_xlabel(_label(x))
    ax.set_ylabel(_label(y))
    if fixed_txt:
        ax.set_title(fixed_txt, fontsize="small")
    return ax


def plot_envelope(envelope, result: QueftsResult | None = None, *, axes=None, **kwargs):
    """Maximum-accumulation / maximum-dilution lines (Janssen 1990 Fig. 1), one panel per nutrient.

    ``envelope`` is the DataFrame of :func:`pyquefts.analysis.uptake_yield_envelope`.
    When ``result`` is given, its (uptake, yield of the modelled organ) pairs are
    overlaid as points; they must fall inside the envelope.
    """
    plt = _plt()
    nuts = list(dict.fromkeys(envelope["nutrient"]))
    if axes is None:
        _, axes = plt.subplots(1, len(nuts), figsize=(4.6 * len(nuts), 3.8), squeeze=False)
        axes = axes[0]
    axes = np.atleast_1d(axes).ravel()
    for a, n in zip(axes, nuts):
        e = envelope[envelope["nutrient"] == n]
        a.plot(e["uptake"], e["Y_acc"], label="max. accumulation", **kwargs)
        a.plot(e["uptake"], e["Y_dil"], label="max. dilution", **kwargs)
        a.fill_between(e["uptake"], e["Y_acc"], e["Y_dil"], alpha=0.12)
        if result is not None:
            upt = _full(result, result[f"{n}_uptake"])
            if "organ" in result:
                organ = _full(result, result["organ"])
                yld = np.full_like(upt, np.nan)
                for code, key in _ORGAN_OUTPUT.items():
                    m = organ == code
                    yld[m] = _full(result, result[key])[m]
            else:
                yld = _full(result, result["store_lim"])
            a.scatter(upt, yld, s=18, color="k", zorder=3, label="model runs")
        a.set_xlabel(f"{n} uptake (kg/ha)")
        a.set_ylabel("yield of modelled organ (kg/ha DM)")
        a.grid(alpha=0.3)
        a.legend(fontsize="small")
    return axes
