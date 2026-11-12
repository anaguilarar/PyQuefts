"""Module 3: spatial QUEFTS simulation over a soil data cube (xarray + dask).

The model core is applied chunk-wise with ``xarray.apply_ufunc``; each chunk is
flattened, run once through the vectorized core and reshaped, so the Python
overhead is per chunk, not per pixel. Memory is bounded by the chunk size.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr

from ..core.model import run_quefts
from ..core.parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from ..core.result import R_OUTPUTS
from ..core.supply import nut_supply1, nut_supply2
from ..crops import get_crop
from .cube import SoilCube

__all__ = ["supply_map", "run_spatial", "response_spatial", "write_netcdf", "write_geotiff"]


def _as_da(value, template: xr.DataArray, name: str) -> xr.DataArray:
    """Broadcast a scalar, 2-D array or DataArray onto the template grid (y, x)."""
    if isinstance(value, xr.DataArray):
        if set(value.dims) != set(template.dims):
            # a field on another grid (e.g. a CHIRTS lat/lon climatology): interpolate onto the cube grid
            from .climate import regrid_to
            lat_dim = next((d for d in value.dims if d.lower() in ("lat", "latitude", "y")), None)
            lon_dim = next((d for d in value.dims if d.lower() in ("lon", "longitude", "x")), None)
            if lat_dim is None or lon_dim is None or value.ndim != 2:
                raise ValueError(f"{name}: cannot regrid a DataArray with dims {value.dims} onto {template.dims}")
            return regrid_to(value, template, lat_dim=lat_dim, lon_dim=lon_dim).rename(name)
        if value.dims != template.dims:
            value = value.transpose(*template.dims)
        if value.shape != template.shape:
            # same dims, coarser grid: nearest-neighbour reindex onto the template grid
            value = value.reindex_like(template, method="nearest")
        return value.rename(name)
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return xr.full_like(template, float(arr), dtype=float).rename(name)
    if arr.shape == template.shape:
        return xr.DataArray(arr, dims=template.dims, coords=template.coords, name=name)
    raise ValueError(f"{name} must be a scalar, a (y, x) array matching the cube subset {template.shape}, or a DataArray")


def supply_map(cube: SoilCube, *, temp=None, supply_method: str = "sattari2014", ptotal_ratio=None,
               bbox=None, mask=None, top: float = 0.0, bottom: float = 20.0, depth_method: str = "weighted",
               clamp_negative: bool | None = None) -> xr.Dataset:
    """Step 1 on a grid: N, P, K base supply maps (kg/ha) plus the model inputs used.

    ``temp`` may be a scalar, a (y, x) array/DataArray, or None if the cube
    maps a ``temp`` variable. ``clamp_negative`` defaults to the behaviour of
    the chosen method (True for Janssen 1990, False for Sattari 2014).
    """
    variables = ["ph", "soc", "kex", "polsen"] + [v for v in ("ptotal", "temp") if v in cube.mapping]
    inp = cube.model_inputs(variables, top, bottom, bbox=bbox, mask=mask, depth_method=depth_method)
    tmpl = inp["ph"]
    ptotal = inp["ptotal"] if "ptotal" in inp else None
    method2 = supply_method.lower().startswith("sattari") or supply_method in ("2", "nutsupply2")
    if method2:
        if temp is None:
            if "temp" not in inp:
                raise ValueError("temp is required for the Sattari (2014) supply method")
            temp_da = inp["temp"]
        else:
            temp_da = _as_da(temp, tmpl, "temp")
        kwargs = dict(ptotal=ptotal, ptotal_ratio=ptotal_ratio)
        func = lambda t, ph, soc, kex, pol, pt=None: nut_supply2(t, ph, soc, kex, pol, ptotal=pt, **{k: v for k, v in kwargs.items() if k != "ptotal"}).to_array()  # noqa: E731
        args = [temp_da, inp["ph"], inp["soc"], inp["kex"], inp["polsen"]] + ([ptotal] if ptotal is not None else [])
    else:
        func = lambda ph, soc, kex, pol, pt=None: nut_supply1(ph, soc, kex, pol, ptotal=pt).to_array()  # noqa: E731
        args = [inp["ph"], inp["soc"], inp["kex"], inp["polsen"]] + ([ptotal] if ptotal is not None else [])

    sup = xr.apply_ufunc(func, *args, output_core_dims=[["nutrient"]], dask="parallelized",
                         output_dtypes=[float], dask_gufunc_kwargs={"output_sizes": {"nutrient": 3}})
    sup = sup.assign_coords(nutrient=["N", "P", "K"])
    out = inp.copy()
    for i, nut in enumerate("NPK"):
        da = sup.isel(nutrient=i, drop=True)
        if clamp_negative is True or (clamp_negative is None and not method2):
            da = da.where(da >= 0, 0.0)
        out[f"{nut}_base_supply"] = da.assign_attrs(units="kg/ha", long_name=f"potential {nut} supply (120-day season)")
    out.attrs.update(supply_method="sattari2014" if method2 else "janssen1990")
    if method2:
        out["temp"] = temp_da.assign_attrs(units="degC", long_name="growing-season mean air temperature")
    return out


def _core_block(Ns, Ps, Ks, Nf, Pf, Kf, leaf, stem, store, season, *, crop: CropParameters, soil: SoilParameters,
                outputs: Sequence[str]) -> np.ndarray:
    """Run the core on flattened blocks; returns (..., n_outputs)."""
    shape = Ns.shape
    sl = soil.replace(N_base_supply=Ns.ravel(), P_base_supply=Ps.ravel(), K_base_supply=Ks.ravel())
    fert = Fertilizer(N=Nf.ravel(), P=Pf.ravel(), K=Kf.ravel())
    biom = Biomass(leaf_att=leaf.ravel(), stem_att=stem.ravel(), store_att=store.ravel(), SeasonLength=season.ravel())
    res = run_quefts(sl, crop, fert, biom, keep_intermediate=any(k not in R_OUTPUTS for k in outputs))
    return np.stack([np.asarray(res[k], dtype=float).reshape(shape) for k in outputs], axis=-1)


def run_spatial(cube: SoilCube | xr.Dataset, crop, *, store_att, leaf_ratio: float, stem_ratio: float,
                fert: Fertilizer | Mapping | None = None, soil: SoilParameters | Mapping | None = None,
                season_length=120.0, temp=None, supply_method: str = "sattari2014", ptotal_ratio=None,
                bbox=None, mask=None, top: float = 0.0, bottom: float = 20.0, depth_method: str = "weighted",
                outputs: Sequence[str] = R_OUTPUTS, keep_inputs: bool = True, compute: bool = False) -> xr.Dataset:
    """Run QUEFTS for every grid cell of a soil cube (or a precomputed supply Dataset).

    Parameters
    ----------
    cube : SoilCube, or a Dataset already holding ``N_base_supply, P_base_supply,
        K_base_supply`` on (y, x) (e.g. from :func:`supply_map`).
    crop : crop name, CropParameters or mapping.
    store_att : attainable storage-organ yield (kg/ha DM): scalar, (y, x) array
        or DataArray (e.g. a water-limited yield map).
    leaf_ratio, stem_ratio : leaf and stem biomass relative to ``store_att``.
    fert : N, P, K applied; scalars or (y, x) maps.
    soil : recoveries and UptakeAdjust (base supplies come from the cube).
    season_length : days; scalar or map.
    temp, supply_method, ptotal_ratio : Step 1 options (see :func:`supply_map`).
    bbox, mask, top, bottom, depth_method : spatial subset and depth aggregation.
    outputs : variables to return (any key of :class:`QueftsResult`).
    keep_inputs : include soil inputs and supplies in the result.
    compute : trigger the dask computation before returning.

    Returns
    -------
    xarray.Dataset on (y, x) with CRS in ``attrs["crs"]`` (and ``rio`` metadata
    when rioxarray is installed).
    """
    if isinstance(crop, str):
        crop = get_crop(crop)
    crop = _coerce(crop, CropParameters, CropParameters)
    soil = _coerce(soil, SoilParameters, SoilParameters)
    fert = _coerce(fert, Fertilizer, Fertilizer)

    if isinstance(cube, SoilCube):
        sup = supply_map(cube, temp=temp, supply_method=supply_method, ptotal_ratio=ptotal_ratio, bbox=bbox,
                         mask=mask, top=top, bottom=bottom, depth_method=depth_method)
        crs = cube.crs
    else:
        sup = cube
        crs = sup.attrs.get("crs", "EPSG:4326")
    tmpl = sup["N_base_supply"]

    store = _as_da(store_att, tmpl, "store_att")
    leaf = (store * leaf_ratio).rename("leaf_att")
    stem = (store * stem_ratio).rename("stem_att")
    Nf, Pf, Kf = (_as_da(getattr(fert, n), tmpl, f"{n}_fert") for n in "NPK")
    season = _as_da(season_length, tmpl, "SeasonLength")

    out = xr.apply_ufunc(
        _core_block, sup["N_base_supply"], sup["P_base_supply"], sup["K_base_supply"], Nf, Pf, Kf, leaf, stem, store,
        season, kwargs={"crop": crop, "soil": soil, "outputs": list(outputs)},
        output_core_dims=[["variable"]], dask="parallelized", output_dtypes=[float],
        dask_gufunc_kwargs={"output_sizes": {"variable": len(outputs)}},
    ).assign_coords(variable=list(outputs))

    res = xr.Dataset({k: out.sel(variable=k, drop=True) for k in outputs})
    units = {"lim": "kg/ha dry matter", "supply": "kg/ha", "uptake": "kg/ha", "gap": "kg/ha"}
    for k in outputs:
        for suf, u in units.items():
            if k.endswith(suf):
                res[k].attrs["units"] = u
    if keep_inputs:
        for k in sup.data_vars:
            res[k] = sup[k]
        res["store_att"] = store.assign_attrs(units="kg/ha dry matter")
    res.attrs.update({
        "crs": str(crs), "crop": crop.name or "default", "model": "QUEFTS (PyQUEFTS port of Rquefts 1.2-8)",
        "supply_method": sup.attrs.get("supply_method", ""), "layer_cm": sup.attrs.get("layer_cm", ""),
        "leaf_ratio": leaf_ratio, "stem_ratio": stem_ratio,
        "fertilizer": f"N={np.asarray(fert.N).mean():g} P={np.asarray(fert.P).mean():g} K={np.asarray(fert.K).mean():g}",
    })
    try:  # georeference for rioxarray users
        import rioxarray  # noqa: F401
        res = res.rio.write_crs(crs, inplace=False)
        res = res.rio.set_spatial_dims(x_dim=tmpl.dims[-1], y_dim=tmpl.dims[-2], inplace=False)
    except Exception:  # pragma: no cover - rioxarray optional / CRS unparsable
        pass
    return res.compute() if compute else res


def response_spatial(cube: SoilCube | xr.Dataset, crop, nutrient: str = "N", rates=None, *,
                     fert: Fertilizer | Mapping | None = None, outputs: Sequence[str] = ("store_lim",),
                     compute: bool = False, **kwargs) -> xr.Dataset:
    """Fertilizer response of one nutrient over a region: ``run_spatial`` stacked along a rate dimension.

    Step 1 (supply maps) is computed once; the model is then run for each
    rate in ``rates`` (other two nutrients held at ``fert``) and the results
    are concatenated along a new dimension named after ``nutrient`` (kg/ha).
    The result is dask-lazy unless ``compute=True``. ``kwargs`` go to
    :func:`run_spatial` (``store_att``, ``leaf_ratio``, ``stem_ratio``,
    ``temp``, ``bbox``, ...).

    Plot with xarray facets, e.g. ``ds.store_lim.plot(col=nutrient)``, or
    with :func:`pyquefts.plot.plot_spatial_response` (area mean and quantile band).
    """
    nutrient = nutrient.upper()
    if nutrient not in ("N", "P", "K"):
        raise ValueError("nutrient must be 'N', 'P' or 'K'")
    rates = np.arange(0.0, 201.0, 10.0) if rates is None else np.asarray(rates, dtype=float).ravel()
    fert = _coerce(fert, Fertilizer, Fertilizer)
    sup_kw = {k: kwargs.pop(k) for k in ("temp", "supply_method", "ptotal_ratio", "bbox", "mask", "top", "bottom",
                                        "depth_method") if k in kwargs}
    if isinstance(cube, SoilCube):
        sup = supply_map(cube, **sup_kw)
        sup.attrs["crs"] = str(cube.crs)
    else:
        sup = cube
    runs = [run_spatial(sup, crop, fert=fert.replace(**{nutrient: float(r)}), outputs=outputs, keep_inputs=False,
                        **kwargs) for r in rates]
    ds = xr.concat(runs, dim=pd.Index(rates, name=nutrient), combine_attrs="drop_conflicts")
    ds[nutrient].attrs.update(units="kg/ha", long_name=f"{nutrient} applied")
    for k in sup.data_vars:
        ds[k] = sup[k]
    ds.attrs.update(runs[0].attrs)
    ds.attrs["fertilizer"] = f"{nutrient} varied; " + ", ".join(
        f"{k}={float(np.mean(getattr(fert, k))):g}" for k in "NPK" if k != nutrient)
    return ds.compute() if compute else ds


def write_netcdf(ds: xr.Dataset, path: str, complevel: int = 4) -> None:
    """Write with zlib compression and float32 storage for all data variables."""
    enc = {k: {"zlib": True, "complevel": complevel, "dtype": "float32"} for k in ds.data_vars}
    ds.to_netcdf(path, encoding=enc)


def write_geotiff(ds: xr.Dataset, path: str, variables: Sequence[str] | None = None, **kwargs) -> None:
    """Write selected variables as a multi-band GeoTIFF (requires rioxarray)."""
    import rioxarray  # noqa: F401

    variables = list(variables or [k for k in ds.data_vars if k in R_OUTPUTS])
    da = ds[variables].to_array("band").astype("float32")
    da = da.rio.write_crs(ds.attrs.get("crs", "EPSG:4326"))
    da = da.rio.write_nodata(np.nan)
    da.attrs["long_name"] = variables
    da.rio.to_raster(path, compress=kwargs.pop("compress", "deflate"), **kwargs)
