"""Parameter containers for the QUEFTS model.

The field names are identical to those used by R-QUEFTS (``quefts_soil()``,
``quefts_crop()``, ``quefts_fert()``, ``quefts_biom()``) so that parameter sets
can be exchanged between the two implementations without renaming.

All numeric fields accept scalars or NumPy-broadcastable arrays. This is what
allows the same core code to run one site, a table of sites or a raster.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields, replace
from typing import Any, Mapping

import numpy as np

ArrayLike = Any  # float | np.ndarray | sequence

#: Default season-length adjustment table of R-QUEFTS: (season length in days,
#: multiplier of the base supply). Base supplies refer to a 120-day season.
UPTAKE_ADJUST_DEFAULT: np.ndarray = np.array(
    [[0.0, 0.0], [40.0, 0.4], [80.0, 0.7], [120.0, 1.0], [240.0, 1.6], [360.0, 2.0], [1000.0, 2.0]]
)


class _ParamBase:
    """Shared helpers for the parameter dataclasses."""

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], strict: bool = False):
        """Build from a mapping; unknown keys are ignored unless ``strict``."""
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if strict and unknown:
            raise KeyError(f"unknown parameters for {cls.__name__}: {sorted(unknown)}")
        return cls(**{k: v for k, v in d.items() if k in names})

    def to_dict(self) -> dict:
        return asdict(self)

    def replace(self, **changes):
        """Return a copy with some fields changed."""
        return replace(self, **changes)

    def __getitem__(self, key: str):
        return getattr(self, key)


@dataclass(frozen=True, eq=False)
class SoilParameters(_ParamBase):
    """Soil nutrient supply and fertilizer recovery (R-QUEFTS ``quefts_soil``).

    Attributes
    ----------
    N_base_supply, P_base_supply, K_base_supply : float or array
        Potential supply of the unfertilized soil (kg/ha) for a standard season
        of 120 days. Estimate with :func:`pyquefts.nut_supply1` or
        :func:`pyquefts.nut_supply2`.
    N_recovery, P_recovery, K_recovery : float or array
        Fraction of applied fertilizer nutrient that the crop can take up.
    UptakeAdjust : ndarray, shape (n, 2)
        Season-length table (days, multiplier), linearly interpolated.
    """

    N_base_supply: ArrayLike = 60.0
    P_base_supply: ArrayLike = 10.0
    K_base_supply: ArrayLike = 60.0
    N_recovery: ArrayLike = 0.5
    P_recovery: ArrayLike = 0.1
    K_recovery: ArrayLike = 0.5
    UptakeAdjust: np.ndarray = field(default_factory=lambda: UPTAKE_ADJUST_DEFAULT.copy())

    def __post_init__(self):
        ua = np.asarray(self.UptakeAdjust, dtype=float)
        if ua.ndim == 1:  # R-style flattened column-major vector (x..., y...)
            ua = ua.reshape(2, -1).T
        if ua.ndim != 2 or ua.shape[1] != 2:
            raise ValueError("UptakeAdjust must be a two-column (days, factor) table")
        if np.any(np.diff(ua[:, 0]) < 0):
            raise ValueError("UptakeAdjust season lengths must be increasing")
        object.__setattr__(self, "UptakeAdjust", ua)


@dataclass(frozen=True, eq=False)
class CropParameters(_ParamBase):
    """Crop nutrient concentration limits (R-QUEFTS ``quefts_crop``).

    ``_minVeg``/``_maxVeg`` and ``_minStore``/``_maxStore`` are the minimum
    (maximum dilution) and maximum (maximum accumulation) mass fractions
    (kg/kg) of N, P and K in vegetative and storage organs. ``Yzero`` is the
    maximum vegetative biomass (kg/ha) at zero storage-organ yield and ``Nfix``
    the fraction of N uptake supplied by biological fixation.

    The defaults are the R-QUEFTS defaults (a maize parameterization).
    """

    NminStore: ArrayLike = 0.0095
    NminVeg: ArrayLike = 0.004
    NmaxStore: ArrayLike = 0.022
    NmaxVeg: ArrayLike = 0.0125
    PminStore: ArrayLike = 0.0017
    PminVeg: ArrayLike = 0.0004
    PmaxStore: ArrayLike = 0.0075
    PmaxVeg: ArrayLike = 0.003
    KminStore: ArrayLike = 0.002
    KminVeg: ArrayLike = 0.005
    KmaxStore: ArrayLike = 0.006
    KmaxVeg: ArrayLike = 0.02
    Yzero: ArrayLike = 400.0
    Nfix: ArrayLike = 0.0
    name: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


@dataclass(frozen=True, eq=False)
class Biomass(_ParamBase):
    """Attainable dry-matter biomass and season length (R-QUEFTS ``quefts_biom``).

    ``leaf_att``, ``stem_att``, ``store_att`` are the attainable (not nutrient
    limited) biomass of leaves, stems and storage organ in kg/ha dry matter;
    ``SeasonLength`` is the growing season in days.
    """

    leaf_att: ArrayLike = 2500.0
    stem_att: ArrayLike = 2500.0
    store_att: ArrayLike = 5000.0
    SeasonLength: ArrayLike = 120.0

    @classmethod
    def from_yield(cls, store_att, leaf_ratio: float, stem_ratio: float, SeasonLength=120.0) -> "Biomass":
        """Build from storage-organ yield and leaf/stem ratios (as in R ``predict``/``batch``)."""
        store_att = np.asarray(store_att, dtype=float)
        return cls(leaf_att=leaf_ratio * store_att, stem_att=stem_ratio * store_att,
                   store_att=store_att, SeasonLength=SeasonLength)


@dataclass(frozen=True, eq=False)
class Fertilizer(_ParamBase):
    """Fertilizer nutrients applied, kg/ha of N, P and K (R-QUEFTS ``quefts_fert``)."""

    N: ArrayLike = 60.0
    P: ArrayLike = 10.0
    K: ArrayLike = 60.0

    @classmethod
    def from_dict(cls, d: Mapping[str, Any], strict: bool = False):
        """Build from a mapping; single-letter keys are case-insensitive (``k`` -> ``K``)."""
        d = {(k.upper() if isinstance(k, str) and len(k) == 1 else k): v for k, v in d.items()}
        return super().from_dict(d, strict=strict)


def _coerce(obj, cls, default_factory):
    """Accept an instance, a mapping or None and return an instance of ``cls``."""
    if obj is None:
        return default_factory()
    if isinstance(obj, cls):
        return obj
    if isinstance(obj, Mapping):
        return cls.from_dict(obj)
    raise TypeError(f"expected {cls.__name__}, mapping or None, got {type(obj).__name__}")
