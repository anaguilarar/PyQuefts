"""Growing-season temperature for the Sattari (2014) N supply from the CHIRTS-ERA5 daily cube.

``nut_supply2`` needs the mean air temperature of the growing season (°C).
This module derives it from the CHIRTS-ERA5 daily Tmax/Tmin Zarr store
published on Hugging Face (Climate Hazards Center product, 0.05°, 1981 to
present; redistributed by A. Aguilar, CC-BY-4.0):

    https://huggingface.co/datasets/aaguilar90/chirts-era5-daily

Only the chunks intersecting the requested region and dates are transferred.
The daily mean is taken as (Tmax + Tmin) / 2, averaged over the season window
and over the requested years (a climatology), then regridded onto the soil
cube grid. No QUEFTS equation is touched; this is an input-preparation helper.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import xarray as xr

CHIRTS_ERA5_URL = "https://huggingface.co/datasets/aaguilar90/chirts-era5-daily/resolve/main/chirts-era5-daily.zarr"

__all__ = ["CHIRTS_ERA5_URL", "open_chirts", "season_mean_temperature", "regrid_to", "sample_at_points",
           "fill_nearest"]


def open_chirts(url: str = CHIRTS_ERA5_URL) -> xr.Dataset:
    """Open the CHIRTS-ERA5 daily Zarr store lazily (requires zarr >= 3 for Zarr v3)."""
    return xr.open_zarr(url, consolidated=True, chunks={})


def _season_slices(start: str, end: str, years: Iterable[int]) -> list[slice]:
    """Time slices for a 'MM-DD' .. 'MM-DD' window in each year (crossing New Year allowed)."""
    out = []
    for y in years:
        s = pd.Timestamp(f"{y}-{start}")
        e = pd.Timestamp(f"{y}-{end}")
        if e < s:
            e = pd.Timestamp(f"{y + 1}-{end}")
        out.append(slice(s, e))
    return out


def season_mean_temperature(ds: xr.Dataset | None = None, *, start: str = "05-15", end: str = "08-31",
                            years: Sequence[int] = tuple(range(1991, 2021)), bbox: tuple | None = None,
                            tmax: str = "Tmax", tmin: str = "Tmin", lat_dim: str = "lat", lon_dim: str = "lon",
                            fill: bool = True) -> xr.DataArray:
    """Mean growing-season air temperature (°C) on the climate grid.

    Parameters
    ----------
    ds : CHIRTS-like dataset with daily ``Tmax``/``Tmin`` (°C) on (time, lat, lon);
        opened from the Hugging Face store when None.
    start, end : season window as ``"MM-DD"``; ``end`` before ``start`` means the
        season crosses New Year.
    years : years whose seasons are averaged (default 1991-2020 climatology).
        Use a single year for a specific season.
    bbox : (xmin, ymin, xmax, ymax) in degrees; strongly recommended for remote stores.
    fill : fill nodata cells (water, outside domain) with the nearest valid
        value so that coastal soil pixels get a temperature.

    Returns
    -------
    DataArray (lat, lon) in °C, attrs describe the window and years.
    """
    if ds is None:
        ds = open_chirts()
    da = ds[[tmax, tmin]]
    if bbox is not None:
        xmin, ymin, xmax, ymax = bbox
        lat = da[lat_dim]
        lat_slice = slice(ymax, ymin) if lat[0] > lat[-1] else slice(ymin, ymax)
        da = da.sel({lat_dim: lat_slice, lon_dim: slice(xmin, xmax)})
        if da.sizes[lat_dim] == 0 or da.sizes[lon_dim] == 0:
            raise ValueError(f"bbox {bbox} does not intersect the climate grid")
    years = list(years)
    parts = [da.sel(time=sl) for sl in _season_slices(start, end, years)]
    parts = [p for p in parts if p.sizes["time"] > 0]
    if not parts:
        raise ValueError(f"no data for {start}..{end} in years {years[0]}..{years[-1]}")
    season = xr.concat(parts, dim="time")
    tmean = ((season[tmax] + season[tmin]) / 2.0).mean("time", skipna=True)
    tmean = tmean.load()
    if fill:
        tmean = fill_nearest(tmean)
    tmean.name = "temp"
    tmean.attrs = {"units": "degC", "long_name": "mean growing-season air temperature",
                   "season": f"{start}..{end}", "years": f"{years[0]}-{years[-1]}", "n_years": len(years),
                   "n_days": int(season.sizes["time"]), "source": ds.attrs.get("title", "CHIRTS-ERA5 daily")}
    return tmean


def fill_nearest(da: xr.DataArray) -> xr.DataArray:
    """Replace NaN cells with the value of the nearest non-NaN cell (2-D arrays)."""
    from scipy import ndimage

    arr = np.asarray(da.values, dtype=float)
    mask = np.isnan(arr)
    if not mask.any() or mask.all():
        return da
    idx = ndimage.distance_transform_edt(mask, return_distances=False, return_indices=True)
    filled = arr[tuple(idx)]
    return da.copy(data=filled)


def regrid_to(da: xr.DataArray, template: xr.DataArray, method: str = "linear", lat_dim: str = "lat",
              lon_dim: str = "lon") -> xr.DataArray:
    """Interpolate a (lat, lon) field onto the (y, x) grid of ``template``.

    ``method`` is ``"linear"`` (bilinear, default) or ``"nearest"``. The result
    carries the template's dims/coords and is chunked like the template when
    the template is dask-backed.
    """
    ydim, xdim = template.dims[-2], template.dims[-1]
    out = da.interp({lat_dim: template[ydim].values, lon_dim: template[xdim].values}, method=method,
                    kwargs={"fill_value": None})
    out = out.rename({lat_dim: ydim, lon_dim: xdim}).assign_coords({ydim: template[ydim], xdim: template[xdim]})
    out = out.transpose(ydim, xdim)
    # edges beyond the climate grid extent: nearest value
    if out.isnull().any():
        out = fill_nearest(out)
    if template.chunks is not None:
        out = out.chunk({ydim: template.chunks[-2][0], xdim: template.chunks[-1][0]})
    out.attrs = dict(da.attrs)
    return out


def sample_at_points(da: xr.DataArray, lon, lat, lat_dim: str = "lat", lon_dim: str = "lon") -> np.ndarray:
    """Values of a (lat, lon) field at coordinates (nearest cell)."""
    lon = xr.DataArray(np.atleast_1d(np.asarray(lon, float)), dims="pt")
    lat = xr.DataArray(np.atleast_1d(np.asarray(lat, float)), dims="pt")
    return da.sel({lat_dim: lat, lon_dim: lon}, method="nearest").values.astype(float)
