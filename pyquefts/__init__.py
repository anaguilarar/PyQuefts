"""PyQUEFTS: a Python implementation of the QUEFTS model (after R-QUEFTS).

QUEFTS (Quantitative Evaluation of the Fertility of Tropical Soils; Janssen et
al., 1990) estimates the native N, P and K supply of a soil from a few topsoil
properties and computes nutrient-limited crop yield from that supply, crop
parameters, fertilizer application and attainable yield.

Three ways to run the model:

* :mod:`pyquefts.point` — user-supplied soil values (``Quefts``, ``run_point``, ``run_batch``)
* :mod:`pyquefts.geo.extract` — coordinates looked up in a soil data cube
* :mod:`pyquefts.geo.spatial` — whole cubes / regions, chunked with dask
"""

from .core import (R_OUTPUTS, mehlich3_to_olsen, p_to_olsen, UPTAKE_ADJUST_DEFAULT, Biomass, CropParameters, Fertilizer, NutrientSupply,
                   QueftsResult, SoilParameters, nut_supply, nut_supply1, nut_supply2, run_quefts)
from .crops import ORGAN_RATIOS, get_crop, list_crops, organ_ratios
from .point import Quefts, run_batch, run_point
from .plot import plot_curves, plot_envelope, plot_response_curves, plot_result, plot_spatial_response, plot_surface

__version__ = "0.1.0"

__all__ = [
    "Quefts", "run_point", "run_batch", "run_quefts",
    "nut_supply", "nut_supply1", "nut_supply2", "NutrientSupply",
    "SoilParameters", "CropParameters", "Biomass", "Fertilizer", "UPTAKE_ADJUST_DEFAULT",
    "QueftsResult", "R_OUTPUTS",
    "get_crop", "list_crops", "organ_ratios", "ORGAN_RATIOS",
    "mehlich3_to_olsen", "p_to_olsen",
    "plot_result", "plot_response_curves", "plot_curves", "plot_spatial_response", "plot_surface", "plot_envelope",
]
