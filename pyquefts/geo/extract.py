"""Module 2: coordinate-based simulation using a soil data cube."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr

from ..core.model import run_quefts
from ..core.parameters import Biomass, Fertilizer, SoilParameters, _coerce
from ..core.result import R_OUTPUTS
from ..core.supply import nut_supply1, nut_supply2
from ..crops import get_crop
from .cube import SoilCube


def _chunk_shape(da: xr.DataArray, y_dim: str, x_dim: str, default: int = 1024) -> tuple[int, int]:
    """Storage chunk size (rows, cols) of a variable: dask chunks, else NetCDF chunking, else a default."""
    if da.chunks is not None:
        return int(da.chunks[da.dims.index(y_dim)][0]), int(da.chunks[da.dims.index(x_dim)][0])
    cs = da.encoding.get("chunksizes")
    if cs:
        return int(cs[da.dims.index(y_dim)]), int(cs[da.dims.index(x_dim)])
    return min(default, da.sizes[y_dim]), min(default, da.sizes[x_dim])


def extract_points(cube: SoilCube, lon: Sequence[float], lat: Sequence[float],
                   variables: Sequence[str] = ("ph", "soc", "kex", "polsen"), top: float = 0.0, bottom: float = 20.0,
                   depth_method: str = "weighted") -> pd.DataFrame:
    """Soil inputs (model units, reference layer) at the grid cells nearest to the coordinates.

    Returns a DataFrame with ``lon, lat, row, col`` and one column per model
    variable; cells outside the cube or with no data are NaN.
    """
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    if lon.shape != lat.shape:
        raise ValueError("lon and lat must have the same length")
    iy, ix = cube.nearest_index(lon, lat)
    valid = (iy >= 0) & (ix >= 0)
    df = pd.DataFrame({"lon": lon, "lat": lat, "row": iy, "col": ix})
    # group the requested cells by storage chunk and read each chunk once
    for k in variables:
        vals = np.full(len(lon), np.nan)
        if valid.any():
            spec = cube.mapping[k]
            da = cube.ds[spec.name]
            cy, cx = _chunk_shape(da, cube.y_dim, cube.x_dim)
            vi = np.nonzero(valid)[0]
            keys = np.stack([iy[vi] // cy, ix[vi] // cx], axis=1)
            for key in np.unique(keys, axis=0):
                sel = vi[(keys[:, 0] == key[0]) & (keys[:, 1] == key[1])]
                r0, c0 = int(key[0] * cy), int(key[1] * cx)
                block = da.isel({cube.y_dim: slice(r0, r0 + cy), cube.x_dim: slice(c0, c0 + cx)})
                block = block.astype("float64") * spec.factor
                if cube.depth_dim and cube.depth_dim in block.dims:
                    block = cube._aggregate_depth(block, top, bottom, depth_method)
                arr = np.asarray(block.values, dtype=float)
                vals[sel] = arr[iy[sel] - r0, ix[sel] - c0]
            if spec.transform is not None:
                vals[vi] = spec.transform(vals[vi])
        df[k] = vals
    return df


def run_points(cube: SoilCube, lon, lat, crop, *, temp=None, store_att=None, leaf_ratio=None, stem_ratio=None,
               biom: Biomass | Mapping | None = None, fert: Fertilizer | Mapping | None = None,
               soil: SoilParameters | Mapping | None = None, supply_method: str = "sattari2014",
               ptotal_ratio=None, top: float = 0.0, bottom: float = 20.0, depth_method: str = "weighted",
               keep_intermediate: bool = False) -> pd.DataFrame:
    """Extract soil properties at coordinates and run QUEFTS at each.

    Parameters
    ----------
    cube : SoilCube
    lon, lat : coordinates (same CRS as the cube)
    crop : crop name, CropParameters or mapping
    temp : growing-season mean temperature: scalar, per-point array, or a
        (lat, lon) DataArray such as :func:`pyquefts.geo.climate.season_mean_temperature`
        (sampled at the nearest climate cell). Required for
        ``supply_method="sattari2014"`` unless the cube maps a ``temp`` variable.
    store_att, leaf_ratio, stem_ratio : attainable storage yield (kg/ha DM,
        scalar or per point) and organ ratios; alternatively give ``biom``.
    fert : N, P, K applied (kg/ha; scalar or per point)
    soil : recoveries / UptakeAdjust (base supplies come from the cube)
    supply_method : ``"sattari2014"`` or ``"janssen1990"``
    ptotal_ratio : total P / P-Olsen ratio when the cube has no total P

    Returns
    -------
    DataFrame: coordinates, soil inputs, supplies and the model outputs.
    """
    if isinstance(crop, str):
        crop = get_crop(crop)
    variables = ["ph", "soc", "kex", "polsen"]
    for opt in ("ptotal", "norg", "temp"):
        if opt in cube.mapping:
            variables.append(opt)
    df = extract_points(cube, lon, lat, variables, top, bottom, depth_method)
    n = len(df)

    if temp is None and "temp" in df:
        temp = df["temp"].to_numpy()
    if isinstance(temp, xr.DataArray):
        from .climate import sample_at_points
        lat_dim = next((d for d in temp.dims if d.lower() in ("lat", "latitude", "y")), "lat")
        lon_dim = next((d for d in temp.dims if d.lower() in ("lon", "longitude", "x")), "lon")
        temp = sample_at_points(temp, lon, lat, lat_dim=lat_dim, lon_dim=lon_dim)
        df["temp"] = temp
    ptotal = df["ptotal"].to_numpy() if "ptotal" in df else None
    if supply_method.lower().startswith("sattari") or supply_method in ("2", "nutsupply2"):
        if temp is None:
            raise ValueError("temp is required for the Sattari (2014) supply method")
        sup = nut_supply2(np.broadcast_to(np.asarray(temp, float), (n,)), df["ph"], df["soc"], df["kex"], df["polsen"],
                          ptotal=ptotal, ptotal_ratio=ptotal_ratio)
    else:
        sup = nut_supply1(df["ph"], df["soc"], df["kex"], df["polsen"], ptotal=ptotal)
    base = _coerce(soil, SoilParameters, SoilParameters)
    soil_p = base.replace(N_base_supply=sup.N_base_supply, P_base_supply=sup.P_base_supply,
                          K_base_supply=sup.K_base_supply)

    if biom is None:
        if store_att is None or leaf_ratio is None or stem_ratio is None:
            raise ValueError("give biom, or store_att with leaf_ratio and stem_ratio")
        biom = Biomass.from_yield(np.broadcast_to(np.asarray(store_att, float), (n,)), leaf_ratio, stem_ratio)
    res = run_quefts(soil_p, crop, fert, biom, keep_intermediate=keep_intermediate)

    for k, v in sup.to_dict().items():
        df[k] = np.broadcast_to(v, (n,))
    for k in res:
        df[k] = np.broadcast_to(res[k], (n,))
    return df


def response_points(cube: SoilCube, lon, lat, crop, nutrient: str = "N", rates=None, *,
                    outputs=("store_lim", "N_uptake", "P_uptake", "K_uptake", "N_gap", "P_gap", "K_gap"),
                    **kwargs) -> pd.DataFrame:
    """Fertilizer response of one nutrient at coordinates looked up in a soil cube.

    The soil is extracted and Step 1 run once per point (:func:`run_points`);
    the model is then run for all points x all ``rates`` in a single
    vectorized call, the other two fertilizer rates held at ``fert``.

    Parameters
    ----------
    nutrient : ``"N"``, ``"P"`` or ``"K"``.
    rates : application rates (kg/ha); default 0..200 step 10.
    outputs : model outputs to return.
    **kwargs : passed to :func:`run_points` (``temp``, ``store_att``,
        ``leaf_ratio``, ``stem_ratio``, ``fert``, ``soil``, ``supply_method``, ...).

    Returns
    -------
    Long DataFrame with one row per (point, rate): ``point`` (0..n-1), ``label``
    ("lon, lat"), ``lon``, ``lat``, the soil inputs and base supplies, the rate
    column named after ``nutrient`` and the requested outputs. Plot it with
    :func:`pyquefts.plot.plot_curves` (``hue="label"``).
    """
    if isinstance(crop, str):
        crop = get_crop(crop)
    nutrient = nutrient.upper()
    if nutrient not in ("N", "P", "K"):
        raise ValueError("nutrient must be 'N', 'P' or 'K'")
    rates = np.arange(0.0, 201.0, 10.0) if rates is None else np.asarray(rates, dtype=float).ravel()
    fert = _coerce(kwargs.pop("fert", None), Fertilizer, Fertilizer)
    soil = _coerce(kwargs.pop("soil", None), SoilParameters, SoilParameters)
    base = run_points(cube, lon, lat, crop, fert=fert, soil=soil, **kwargs)
    n, m = len(base), rates.size

    col = lambda name: base[name].to_numpy(dtype=float)[:, None]  # noqa: E731  (n, 1)
    soil_p = soil.replace(N_base_supply=col("N_base_supply"), P_base_supply=col("P_base_supply"),
                          K_base_supply=col("K_base_supply"))
    fert_vals = {k: np.broadcast_to(np.asarray(getattr(fert, k), float).reshape(-1, 1), (n, 1)) for k in "NPK"}
    fert_vals[nutrient] = rates[None, :]
    # same biomass as run_points built, shaped (n, 1)
    if kwargs.get("biom") is None:
        store = np.broadcast_to(np.asarray(kwargs["store_att"], float), (n,))[:, None]
        biom_p = Biomass.from_yield(store, kwargs["leaf_ratio"], kwargs["stem_ratio"])
    else:
        b = _coerce(kwargs["biom"], Biomass, Biomass)
        biom_p = Biomass(**{k: np.broadcast_to(np.asarray(getattr(b, k), float).reshape(-1, 1), (n, 1))
                            for k in ("leaf_att", "stem_att", "store_att", "SeasonLength")})
    res = run_quefts(soil_p, crop, Fertilizer(**fert_vals), biom_p, keep_intermediate=any(k not in R_OUTPUTS for k in outputs))

    keep = [c for c in base.columns if c not in R_OUTPUTS]
    out = base[keep].loc[np.repeat(np.arange(n), m)].reset_index(drop=True)
    out.insert(0, "label", [f"{x:.3f}, {y:.3f}" for x, y in zip(out["lon"], out["lat"])])
    out.insert(0, "point", np.repeat(np.arange(n), m))
    out["store_att"] = np.broadcast_to(biom_p.store_att, (n, m)).ravel()
    out[nutrient] = np.tile(rates, n)
    for k in "NPK":
        if k != nutrient:
            out[f"{k}_applied"] = np.broadcast_to(fert_vals[k], (n, m)).ravel()
    for k in outputs:
        out[k] = np.broadcast_to(np.asarray(res[k], float), (n, m)).ravel()
    out.attrs.update(nutrient=nutrient, crop=crop.name or "default")
    return out
