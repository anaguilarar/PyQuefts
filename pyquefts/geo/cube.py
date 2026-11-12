"""Soil data cubes: opening, variable mapping, unit conversion and depth aggregation.

A :class:`SoilCube` wraps an ``xarray.Dataset`` (NetCDF, Zarr or a stack of
GeoTIFFs) and knows how to produce the five QUEFTS soil inputs in model units
(pH, SOC g/kg, Kex mmol/kg, P-Olsen mg/kg, optionally Norg g/kg and total P
mg/kg) for a reference soil layer (default 0-20 cm, the layer QUEFTS was
calibrated on).

Nothing here assumes a naming convention: the mapping between dataset
variables and model variables is explicit (:class:`VariableSpec`) and is
inspected from the dataset metadata where possible.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional, Sequence

import numpy as np
import xarray as xr

from ..core.conversions import P_METHODS, p_to_olsen

#: Model variables that the cube may provide and their model units.
MODEL_VARS = {
    "ph": "pH (H2O)",
    "soc": "g/kg organic carbon",
    "kex": "mmol/kg exchangeable K",
    "polsen": "mg/kg P-Olsen",
    "norg": "g/kg organic N (optional)",
    "ptotal": "mg/kg total P (optional)",
    "temp": "degrees C (optional, growing-season mean)",
}

#: Scale factors from common storage units to model units.
UNIT_FACTORS = {
    # pH
    ("ph", "ph*10"): 0.1, ("ph", "ph"): 1.0, ("ph", "-"): 1.0, ("ph", ""): 1.0,
    # organic carbon -> g/kg
    ("soc", "dg/kg"): 0.1, ("soc", "g/kg"): 1.0, ("soc", "%"): 10.0, ("soc", "cg/kg"): 0.01, ("soc", "g/100g"): 10.0,
    # exchangeable K -> mmol/kg  (1 cmol(+)/kg = 10 mmol(+)/kg)
    ("kex", "cmol(+)/kg"): 10.0, ("kex", "cmol/kg"): 10.0, ("kex", "cmolc/kg"): 10.0, ("kex", "mmol/kg"): 1.0,
    ("kex", "mmol(+)/kg"): 1.0, ("kex", "mg/kg"): 1.0 / 39.0983,  # mg K/kg -> mmol/kg
    # P -> mg/kg
    ("polsen", "mg/kg"): 1.0, ("polsen", "ppm"): 1.0, ("ptotal", "mg/kg"): 1.0, ("ptotal", "ppm"): 1.0,
    # organic N -> g/kg
    ("norg", "cg/kg"): 0.01, ("norg", "g/kg"): 1.0, ("norg", "%"): 10.0, ("norg", "mg/kg"): 0.001,
    # temperature
    ("temp", "degc"): 1.0, ("temp", "c"): 1.0, ("temp", "degrees_celsius"): 1.0, ("temp", "celsius"): 1.0,
}

#: Candidate dataset variable names for each model variable (lower-case).
NAME_CANDIDATES = {
    "ph": ["phh2o", "ph_h2o", "ph", "phw", "ph_water", "soil_ph"],
    "soc": ["soc", "orgc", "oc", "org_c", "organic_carbon", "carbon"],
    "kex": ["k", "kex", "k_exch", "exch_k", "exchangeable_k", "kexch", "potassium"],
    "polsen": ["p", "polsen", "p_olsen", "olsen_p", "pex", "p_avail", "phosphorus"],
    "norg": ["nitrogen", "n", "norg", "ntot", "n_tot", "total_n"],
    "ptotal": ["ptot", "ptotal", "p_total", "total_p"],
    "temp": ["temp", "tavg", "tmean", "temperature"],
}


@dataclass
class VariableSpec:
    """How one model variable is obtained from the dataset.

    ``name`` is the dataset variable; ``factor`` multiplies the stored values to
    get model units; ``transform`` (optional) is applied after the factor and
    may be used, for example, to convert Mehlich-3 P to P-Olsen.
    """

    name: str
    factor: float = 1.0
    transform: Optional[Callable[[np.ndarray], np.ndarray]] = None
    source_units: str = ""


def _norm_unit(u: str) -> str:
    return re.sub(r"\s+", "", str(u)).lower().replace("°", "deg")


def _parse_depth_label(label) -> tuple[float, float]:
    """'0-5', '5-15cm', '0_5' ... -> (top, bottom) in cm. Numeric labels are treated as midpoints."""
    s = str(label).lower().replace("cm", "").replace("_", "-").strip()
    m = re.match(r"^\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*$", s)
    if m:
        return float(m.group(1)), float(m.group(2))
    try:
        v = float(s)
        return v, v
    except ValueError as e:
        raise ValueError(f"cannot parse depth label {label!r}") from e


def depth_weights(depth_labels: Sequence, top: float = 0.0, bottom: float = 20.0) -> np.ndarray:
    """Thickness weights of each layer inside the target interval ``[top, bottom]`` cm.

    Layers are given as labels like ``'0-5', '5-15', '15-30'``. The weights are
    the overlap thickness with the target interval, normalized to sum to one.
    """
    w = []
    for lab in depth_labels:
        t, b = _parse_depth_label(lab)
        if t == b:  # point depth: unit weight if inside
            w.append(1.0 if top <= t <= bottom else 0.0)
        else:
            w.append(max(0.0, min(b, bottom) - max(t, top)))
    w = np.asarray(w, dtype=float)
    if w.sum() <= 0:
        raise ValueError(f"no depth layer overlaps {top}-{bottom} cm: {list(depth_labels)}")
    return w / w.sum()


class SoilCube:
    """A soil data cube with explicit mapping to QUEFTS inputs.

    Parameters
    ----------
    ds : xarray.Dataset (lazily loaded; use ``chunks`` when opening).
    mapping : dict model_var -> VariableSpec | str. Strings are dataset
        variable names with the unit factor inferred from the ``units``
        attribute. Missing entries are auto-detected from
        :data:`NAME_CANDIDATES`.
    depth_dim, x_dim, y_dim : dimension names (auto-detected if None).
    crs : CRS string/EPSG for output georeferencing (read from attrs if None).
    p_method : soil-test method of the P variable: ``"olsen"`` (model input,
        default) or ``"mehlich3"`` (converted to Olsen P with McDowell et al.
        2023, see :mod:`pyquefts.core.conversions`). Ignored when ``mapping``
        gives an explicit ``VariableSpec`` for ``polsen``.
    calcareous : use the calcareous-soil P conversion (default False).
    """

    def __init__(self, ds: xr.Dataset, mapping: Mapping[str, VariableSpec | str] | None = None,
                 depth_dim: str | None = None, x_dim: str | None = None, y_dim: str | None = None,
                 crs: str | None = None, p_method: str = "olsen", calcareous: bool = False):
        self.ds = ds
        self.p_method = p_method.lower().replace("-", "").replace(" ", "")
        if self.p_method not in P_METHODS:
            raise ValueError(f"unknown p_method {p_method!r}; choose from {sorted(P_METHODS)}")
        self.calcareous = calcareous
        self.x_dim = x_dim or self._guess_dim(("x", "lon", "longitude"))
        self.y_dim = y_dim or self._guess_dim(("y", "lat", "latitude"))
        self.depth_dim = depth_dim or self._guess_dim(("depth", "depths", "z", "layer"), required=False)
        self.crs = crs or self._guess_crs()
        self.mapping: dict[str, VariableSpec] = {}
        for k, v in (mapping or {}).items():
            self.mapping[k] = v if isinstance(v, VariableSpec) else self._spec_from_name(k, v)
        for k in MODEL_VARS:
            if k not in self.mapping:
                spec = self._auto_spec(k)
                if spec is not None:
                    self.mapping[k] = spec
        if self.p_method != "olsen" and "polsen" in self.mapping and self.mapping["polsen"].transform is None:
            spec = self.mapping["polsen"]
            conv = (lambda a, _m=self.p_method, _c=self.calcareous: p_to_olsen(a, _m, calcareous=_c))
            self.mapping["polsen"] = VariableSpec(spec.name, spec.factor, conv, spec.source_units + f" ({self.p_method})")

    # ------------------------------------------------------------ constructors
    @classmethod
    def open(cls, path: str, chunks: Mapping | str | None = "auto", engine: str | None = None, **kwargs) -> "SoilCube":
        """Open NetCDF (``.nc``) or Zarr (``.zarr``) lazily."""
        p = str(path)
        if p.endswith(".zarr") or p.rstrip("/").endswith(".zarr"):
            ds = xr.open_zarr(p, chunks=chunks if chunks != "auto" else None)
        else:
            ds = xr.open_dataset(p, chunks=chunks, engine=engine, decode_times=False)
        return cls(ds, **kwargs)

    @classmethod
    def from_rasters(cls, files: Mapping[str, str], **kwargs) -> "SoilCube":
        """Build a single-layer cube from GeoTIFFs: ``{model_var: path}`` (requires rioxarray)."""
        import rioxarray  # noqa: F401

        arrays = {}
        for k, f in files.items():
            da = rioxarray.open_rasterio(f, masked=True).squeeze("band", drop=True)
            arrays[k] = da
        ds = xr.Dataset(arrays)
        mapping = {k: VariableSpec(k, 1.0) for k in files}
        crs = kwargs.pop("crs", None) or str(next(iter(arrays.values())).rio.crs)
        return cls(ds, mapping=mapping, x_dim="x", y_dim="y", crs=crs, **kwargs)

    # ------------------------------------------------------------ helpers
    def _guess_dim(self, candidates, required=True):
        dims = {d.lower(): d for d in self.ds.dims}
        for c in candidates:
            if c in dims:
                return dims[c]
        if required:
            raise ValueError(f"cannot find a dimension among {candidates} in {list(self.ds.dims)}")
        return None

    def _guess_crs(self):
        for key in ("crs", "spatial_ref", "crs_wkt", "proj4"):
            if key in self.ds.attrs:
                return self.ds.attrs[key]
        if "spatial_ref" in self.ds.coords:
            return self.ds["spatial_ref"].attrs.get("crs_wkt", "EPSG:4326")
        return "EPSG:4326"

    def _spec_from_name(self, model_var: str, name: str) -> VariableSpec:
        if name not in self.ds:
            raise KeyError(f"variable {name!r} not in dataset")
        units = _norm_unit(self.ds[name].attrs.get("units", ""))
        factor = UNIT_FACTORS.get((model_var, units))
        if factor is None:
            raise ValueError(f"unknown units {units!r} for {model_var} (variable {name!r}); "
                             f"pass a VariableSpec with an explicit factor")
        return VariableSpec(name, factor, source_units=units)

    def _auto_spec(self, model_var: str) -> VariableSpec | None:
        names = {v.lower(): v for v in self.ds.data_vars}
        for cand in NAME_CANDIDATES[model_var]:
            if cand in names:
                try:
                    return self._spec_from_name(model_var, names[cand])
                except ValueError:
                    return None
        return None

    def describe(self) -> str:
        lines = [f"SoilCube dims: x={self.x_dim} ({self.ds.sizes[self.x_dim]}), y={self.y_dim} "
                 f"({self.ds.sizes[self.y_dim]}), depth={self.depth_dim} "
                 f"({list(self.ds[self.depth_dim].values) if self.depth_dim else 'none'})", f"crs: {str(self.crs)[:60]}..."]
        for k, spec in self.mapping.items():
            lines.append(f"  {k:8s} <- {spec.name!r} [{spec.source_units}] x {spec.factor:g}"
                         + (" + transform" if spec.transform else ""))
        if self.p_method != "olsen":
            lines.append(f"  P method: {self.p_method} -> Olsen (McDowell et al. 2023, "
                         f"{'calcareous' if self.calcareous else 'non-calcareous'} equation)")
        missing = [k for k in ("ph", "soc", "kex", "polsen") if k not in self.mapping]
        if missing:
            lines.append(f"  MISSING required: {missing}")
        return "\n".join(lines)

    def __repr__(self):
        return self.describe()

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        x, y = self.ds[self.x_dim].values, self.ds[self.y_dim].values
        return float(x.min()), float(y.min()), float(x.max()), float(y.max())

    # ------------------------------------------------------------ extraction
    def model_inputs(self, variables: Sequence[str] = ("ph", "soc", "kex", "polsen"), top: float = 0.0,
                     bottom: float = 20.0, bbox: tuple | None = None, mask: xr.DataArray | None = None,
                     depth_method: str = "weighted") -> xr.Dataset:
        """Model-unit inputs for the reference layer as a 2-D (y, x) Dataset.

        Parameters
        ----------
        variables : which model variables to produce.
        top, bottom : target soil layer (cm) for depth aggregation.
        bbox : (xmin, ymin, xmax, ymax) subset in coordinate units.
        mask : boolean DataArray on (y, x); False cells become NaN.
        depth_method : ``"weighted"`` thickness-weighted mean over the layers
            overlapping [top, bottom]; ``"top"`` first layer only; ``"mean"``
            plain mean of overlapping layers.
        """
        ds = self.ds
        if bbox is not None:
            ds = self._subset_bbox(ds, bbox)
        out = {}
        for k in variables:
            if k not in self.mapping:
                raise KeyError(f"model variable {k!r} is not mapped in this cube: {self.describe()}")
            spec = self.mapping[k]
            da = ds[spec.name].astype("float64") * spec.factor
            if self.depth_dim and self.depth_dim in da.dims:
                da = self._aggregate_depth(da, top, bottom, depth_method)
            if spec.transform is not None:
                da = xr.apply_ufunc(spec.transform, da, dask="parallelized", output_dtypes=[float])
            da.attrs = {"units": MODEL_VARS[k], "source": spec.name}
            out[k] = da
        res = xr.Dataset(out)
        if mask is not None:
            res = res.where(mask)
        res.attrs["layer_cm"] = f"{top}-{bottom}"
        res.attrs["crs"] = str(self.crs)
        res.attrs["p_method"] = self.p_method
        return res

    def _subset_bbox(self, ds, bbox):
        xmin, ymin, xmax, ymax = bbox
        x, y = ds[self.x_dim], ds[self.y_dim]
        xs = slice(xmin, xmax) if x[0] < x[-1] else slice(xmax, xmin)
        ys = slice(ymin, ymax) if y[0] < y[-1] else slice(ymax, ymin)
        sub = ds.sel({self.x_dim: xs, self.y_dim: ys})
        if sub.sizes[self.x_dim] == 0 or sub.sizes[self.y_dim] == 0:
            raise ValueError(f"bbox {bbox} does not intersect the cube bounds {self.bounds}")
        return sub

    def _aggregate_depth(self, da: xr.DataArray, top: float, bottom: float, method: str) -> xr.DataArray:
        labels = list(da[self.depth_dim].values)
        w = depth_weights(labels, top, bottom)
        sel = np.nonzero(w)[0]
        if method == "top":
            return da.isel({self.depth_dim: int(sel[0])}, drop=True)
        sub = da.isel({self.depth_dim: sel})
        if method == "mean":
            return sub.mean(self.depth_dim, skipna=False)
        if method != "weighted":
            raise ValueError("depth_method must be 'weighted', 'top' or 'mean'")
        wda = xr.DataArray(w[sel], dims=(self.depth_dim,), coords={self.depth_dim: sub[self.depth_dim]})
        return (sub * wda).sum(self.depth_dim, skipna=False)

    def nearest_index(self, lon, lat) -> tuple[np.ndarray, np.ndarray]:
        """Row/column indices of the grid cells nearest to the coordinates (NaN-safe, no interpolation)."""
        x = self.ds[self.x_dim].values
        y = self.ds[self.y_dim].values
        lon = np.atleast_1d(np.asarray(lon, dtype=float))
        lat = np.atleast_1d(np.asarray(lat, dtype=float))
        dx = np.abs(np.median(np.diff(x)))
        dy = np.abs(np.median(np.diff(y)))
        ix = np.rint((lon - x[0]) / (x[1] - x[0])).astype(int)
        iy = np.rint((lat - y[0]) / (y[1] - y[0])).astype(int)
        inside = (lon >= x.min() - dx / 2) & (lon <= x.max() + dx / 2) & (lat >= y.min() - dy / 2) & (lat <= y.max() + dy / 2)
        ix = np.where(inside, np.clip(ix, 0, len(x) - 1), -1)
        iy = np.where(inside, np.clip(iy, 0, len(y) - 1), -1)
        return iy, ix
