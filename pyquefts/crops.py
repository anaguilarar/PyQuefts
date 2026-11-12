"""Crop parameter sets shipped with R-QUEFTS (``inst/extdata/quefts_crop_pars.csv``)."""

from __future__ import annotations

import csv
from functools import lru_cache
from importlib import resources

from .core.parameters import CropParameters

_NUMERIC = ("NminStore", "NminVeg", "NmaxStore", "NmaxVeg", "PminStore", "PminVeg", "PmaxStore", "PmaxVeg",
            "KminStore", "KminVeg", "KmaxStore", "KmaxVeg", "Yzero", "Nfix")

#: Leaf and stem dry matter relative to the storage organ, as documented in the
#: Rquefts manual (``predict``/``batch`` examples). Other crops need user values.
ORGAN_RATIOS = {
    "maize": (0.46, 0.56),
    "potato": (0.17, 0.14),
    "rice": (0.18, 0.67),
}


@lru_cache(maxsize=1)
def _table() -> dict[str, dict]:
    with resources.files("pyquefts.data").joinpath("quefts_crop_pars.csv").open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        d = {k: float(r[k]) for k in _NUMERIC}
        d["name"] = r["crop"]
        out[r["crop"].lower()] = d
    return out


def list_crops() -> list[str]:
    """Names of the available crop parameter sets."""
    return [v["name"] for v in _table().values()]


def get_crop(name: str = "") -> CropParameters:
    """Crop parameters by name (case-insensitive), as R ``quefts_crop(name)``.

    An empty name returns the R-QUEFTS default parameter set (maize values).
    """
    if not name:
        return CropParameters()
    key = name.strip().lower()
    tbl = _table()
    if key not in tbl:
        raise KeyError(f"crop {name!r} not found. Choose from: {', '.join(list_crops())}")
    return CropParameters(**tbl[key])


def organ_ratios(name: str) -> tuple[float, float]:
    """Documented (leaf, stem) ratios relative to storage-organ yield for a few crops."""
    key = name.strip().lower()
    if key not in ORGAN_RATIOS:
        raise KeyError(f"no documented organ ratios for {name!r}; provide leaf_ratio and stem_ratio explicitly")
    return ORGAN_RATIOS[key]
