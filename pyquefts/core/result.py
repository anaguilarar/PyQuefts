"""Container for QUEFTS outputs."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Iterator

import numpy as np

#: The 12 outputs returned by R-QUEFTS ``run()``, in R order.
R_OUTPUTS = ("leaf_lim", "stem_lim", "store_lim", "N_supply", "P_supply", "K_supply",
             "N_uptake", "P_uptake", "K_uptake", "N_gap", "P_gap", "K_gap")


class QueftsResult(Mapping):
    """Dict-like, read-only container of model outputs (NumPy arrays of a common shape).

    ``res["store_lim"]`` returns the array; for scalar runs use ``res.scalar("store_lim")``
    or ``res.to_dict()`` which converts 0-d arrays to floats.
    """

    def __init__(self, data: dict, inputs: Mapping | None = None):
        self._data = dict(data)
        shapes = {np.shape(v) for v in self._data.values()}
        self.shape = max(shapes, key=len) if shapes else ()
        #: inputs of the run (fertilizer, biomass, base supplies, recoveries), recorded by
        #: :meth:`pyquefts.Quefts.run` so that :meth:`plot` can find the one that varies.
        self.inputs: dict | None = dict(inputs) if inputs is not None else None

    # Mapping interface -----------------------------------------------------
    def __getitem__(self, key: str) -> np.ndarray:
        return self._data[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __repr__(self) -> str:
        if self.shape == ():
            body = ", ".join(f"{k}={float(self._data[k]):.4g}" for k in R_OUTPUTS if k in self._data)
        else:
            body = f"shape={self.shape}, {len(self._data)} variables"
        return f"QueftsResult({body})"

    # Conversions -----------------------------------------------------------
    def scalar(self, key: str) -> float:
        return float(np.asarray(self._data[key]).ravel()[0])

    def to_dict(self, r_outputs_only: bool = False) -> dict:
        keys = R_OUTPUTS if r_outputs_only else self._data.keys()
        out = {}
        for k in keys:
            v = np.asarray(self._data[k])
            out[k] = v.item() if v.ndim == 0 else v
        return out

    def to_array(self, keys=R_OUTPUTS) -> np.ndarray:
        """Stack selected outputs along the last axis."""
        return np.stack([np.asarray(self._data[k], dtype=float) for k in keys], axis=-1)

    def to_dataframe(self, r_outputs_only: bool = False, index=None, include_inputs: bool = False):
        """Flatten to a pandas DataFrame (one row per element).

        ``include_inputs=True`` prepends the recorded inputs that vary across
        elements (e.g. the N rates of a response curve).
        """
        import pandas as pd
        keys = R_OUTPUTS if r_outputs_only else list(self._data.keys())
        cols = {}
        if include_inputs and self.inputs:
            for k, v in self.inputs.items():
                if np.ndim(v) > 0 and np.size(v) > 1:
                    cols[k] = np.broadcast_to(np.asarray(v, dtype=float), self.shape).ravel()
        cols.update({k: np.asarray(self._data[k]).ravel() for k in keys})
        return pd.DataFrame(cols, index=index)

    # Plotting (extension, requires matplotlib) ------------------------------
    def plot(self, x: str | None = None, y="store_lim", *, hue: str | None = None, ax=None, **kwargs):
        """Plot outputs against the input that varies; see :func:`pyquefts.plot.plot_result`.

        >>> q = Quefts(fert=dict(N=np.arange(0, 160, 10), P=10, K=60))
        >>> ax = q.run().plot()                        # store_lim vs N
        >>> axes = q.run().plot(y=["store_lim", "N_uptake", "N_gap"])
        """
        from ..plot import plot_result
        return plot_result(self, x=x, y=y, hue=hue, ax=ax, **kwargs)

    def to_dataset(self, coords: Mapping, dims: tuple, attrs: Mapping | None = None, r_outputs_only=False):
        """Wrap outputs as an ``xarray.Dataset`` on the given dims/coords."""
        import xarray as xr
        keys = R_OUTPUTS if r_outputs_only else list(self._data.keys())
        ds = xr.Dataset({k: (dims, np.asarray(self._data[k])) for k in keys}, coords=dict(coords))
        if attrs:
            ds.attrs.update(attrs)
        return ds
