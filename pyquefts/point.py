"""Module 1: point-based simulation with user-supplied values.

Two entry points:

* :class:`Quefts` — an object with the same shape as the R-QUEFTS model
  (``soil``, ``crop``, ``fert``, ``biom`` and ``run()``), convenient for
  interactive use and parameter sweeps.
* :func:`run_batch` — the vectorized equivalent of R ``batch()``: many
  supply/fertilizer/attainable-yield combinations in one call.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .core.model import run_quefts
from .core.parameters import Biomass, CropParameters, Fertilizer, SoilParameters, _coerce
from .core.result import QueftsResult
from .crops import get_crop


class Quefts:
    """A QUEFTS model instance (R-QUEFTS ``quefts()`` equivalent).

    Parameters
    ----------
    soil : SoilParameters or mapping (default: R-QUEFTS example soil)
    crop : CropParameters, mapping or crop name (default: maize-like defaults)
    fert : Fertilizer or mapping with N, P, K (default: 60, 10, 60 kg/ha)
    biom : Biomass or mapping (default: 2500/2500/5000 kg/ha, 120 days)

    Examples
    --------
    >>> q = Quefts(crop="barley", fert={"N": 100, "P": 0, "K": 0},
    ...            biom={"leaf_att": 2200, "stem_att": 2700, "store_att": 4800, "SeasonLength": 110})
    >>> round(q.run().scalar("store_lim"), 1)
    2877.4
    """

    def __init__(self, soil=None, crop=None, fert=None, biom=None):
        if isinstance(crop, str):
            crop = get_crop(crop)
        self.soil: SoilParameters = _coerce(soil, SoilParameters, SoilParameters)
        self.crop: CropParameters = _coerce(crop, CropParameters, CropParameters)
        self.fert: Fertilizer = _coerce(fert, Fertilizer, Fertilizer)
        self.biom: Biomass = _coerce(biom, Biomass, Biomass)

    # ---- convenience setters (R: q$N <- 100, q["soil$N_recovery"] <- .6)
    def set(self, **kwargs) -> "Quefts":
        """Update fields on any of the four parameter blocks, e.g. ``q.set(N=100, SeasonLength=110)``."""
        for k, v in kwargs.items():
            for block in ("soil", "crop", "fert", "biom"):
                obj = getattr(self, block)
                if hasattr(obj, k):
                    setattr(self, block, obj.replace(**{k: v}))
                    break
            else:
                raise KeyError(f"unknown parameter {k!r}")
        return self

    def get(self, key: str):
        for block in ("soil", "crop", "fert", "biom"):
            obj = getattr(self, block)
            if hasattr(obj, key):
                return getattr(obj, key)
        raise KeyError(key)

    def __getitem__(self, key: str):
        return self.get(key)

    def run(self, keep_intermediate: bool = True) -> QueftsResult:
        """Run the model with the current parameters.

        Array-valued parameters (e.g. ``fert=dict(N=[0, 50, 100], P=10, K=60)``)
        give array outputs, the vectorized equivalent of the R loop
        ``for (i in 1:length(N)) { q["N"] <- N[i]; run(q) }``. The inputs are
        recorded on the result so that ``result.plot()`` knows the x axis.
        """
        res = run_quefts(self.soil, self.crop, self.fert, self.biom, keep_intermediate=keep_intermediate)
        res.inputs = self.inputs()
        return res

    def inputs(self) -> dict:
        """The varying-capable inputs of this instance (fertilizer, biomass, base supplies, recoveries)."""
        out = {k: getattr(self.fert, k) for k in ("N", "P", "K")}
        out.update({k: getattr(self.biom, k) for k in ("leaf_att", "stem_att", "store_att", "SeasonLength")})
        out.update({k: getattr(self.soil, k) for k in ("N_base_supply", "P_base_supply", "K_base_supply",
                                                        "N_recovery", "P_recovery", "K_recovery")})
        return out

    def response(self, nutrient: str, rates, outputs=("store_lim", "N_uptake", "P_uptake", "K_uptake")):
        """Response table of one nutrient (the other two held at ``self.fert``); see :func:`pyquefts.analysis.response_curve`."""
        from .analysis import response_curve
        return response_curve(nutrient, rates, self.crop, soil=self.soil, fert=self.fert, biom=self.biom, outputs=outputs)

    def plot(self, rates=None, *, nutrients="NPK", y: str = "store_lim", axes=None, **kwargs):
        """N, P and K response panels for this instance; see :func:`pyquefts.plot.plot_response_curves`."""
        from .plot import plot_response_curves
        return plot_response_curves(self, rates, nutrients=nutrients, y=y, axes=axes, **kwargs)

    def __repr__(self) -> str:
        return (f"Quefts(crop={self.crop.name or 'default'!r}, fert=({self.fert.N}, {self.fert.P}, {self.fert.K}), "
                f"biom=({self.biom.leaf_att}, {self.biom.stem_att}, {self.biom.store_att}, {self.biom.SeasonLength}))")


def run_batch(crop, supply, fert, yatt, leaf_ratio: float, stem_ratio: float, *, soil=None,
              season_length=120.0, var: str = "all") -> QueftsResult | np.ndarray:
    """Vectorized batch run (R-QUEFTS ``batch()``).

    Parameters
    ----------
    crop : CropParameters, mapping or crop name
    supply : array-like (n, 3) of base N, P, K supply (kg/ha); rows broadcast
    fert : array-like (n, 3) of N, P, K fertilizer (kg/ha); rows broadcast
    yatt : attainable storage-organ yield (kg/ha), shape (n,) or scalar
    leaf_ratio, stem_ratio : leaf and stem biomass relative to ``yatt``
    soil : SoilParameters for recoveries and UptakeAdjust (base supplies are
        overridden by ``supply``)
    season_length : days
    var : ``"yield"`` returns ``store_lim``; ``"gap"`` returns (n, 3) gaps;
        ``"all"`` returns the full :class:`QueftsResult`.
    """
    if isinstance(crop, str):
        crop = get_crop(crop)
    base = _coerce(soil, SoilParameters, SoilParameters)
    supply = np.atleast_2d(np.asarray(supply, dtype=float))
    fert = np.atleast_2d(np.asarray(fert, dtype=float))
    yatt = np.asarray(yatt, dtype=float)
    n = max(supply.shape[0], fert.shape[0], yatt.size)
    supply = np.broadcast_to(supply, (n, 3))
    fert = np.broadcast_to(fert, (n, 3))
    yatt = np.broadcast_to(yatt.ravel() if yatt.ndim else yatt, (n,))

    soil_p = base.replace(N_base_supply=supply[:, 0], P_base_supply=supply[:, 1], K_base_supply=supply[:, 2])
    fert_p = Fertilizer(N=fert[:, 0], P=fert[:, 1], K=fert[:, 2])
    biom_p = Biomass.from_yield(yatt, leaf_ratio, stem_ratio, SeasonLength=season_length)
    res = run_quefts(soil_p, crop, fert_p, biom_p, keep_intermediate=(var == "all"))
    if var == "yield":
        return res["store_lim"]
    if var == "gap":
        return res.to_array(("N_gap", "P_gap", "K_gap"))
    if var == "all":
        return res
    raise ValueError("var must be 'yield', 'gap' or 'all'")


def run_point(ph=None, soc=None, kex=None, polsen=None, ptotal=None, temp=None, *, crop="", fert=None,
              biom=None, supply_method: str = "sattari2014", ptotal_ratio=None, recovery: Mapping | None = None,
              uptake_adjust=None, p_method: str = "olsen", calcareous: bool = False) -> QueftsResult:
    """One-call point simulation from raw soil properties.

    Runs Step 1 (:func:`nut_supply1` or :func:`nut_supply2`) and then the
    model. ``temp`` is required for ``supply_method="sattari2014"``.
    ``polsen`` is taken as Olsen P unless ``p_method="mehlich3"``, in which
    case it is converted first (McDowell et al. 2023).
    """
    from .core.conversions import p_to_olsen
    from .core.supply import nut_supply1, nut_supply2

    polsen = p_to_olsen(polsen, p_method, calcareous=calcareous)

    if supply_method.lower() in ("sattari2014", "2", "sattari", "nutsupply2"):
        if temp is None:
            raise ValueError("temp (growing-season mean temperature, C) is required for the Sattari 2014 supply")
        sup = nut_supply2(temp, ph, soc, kex, polsen, ptotal=ptotal, ptotal_ratio=ptotal_ratio)
    else:
        sup = nut_supply1(ph, soc, kex, polsen, ptotal=ptotal)
    soil_kwargs = dict(recovery or {})
    if uptake_adjust is not None:
        soil_kwargs["UptakeAdjust"] = uptake_adjust
    soil = sup.to_soil(**soil_kwargs)
    return Quefts(soil=soil, crop=crop, fert=fert, biom=biom).run()
