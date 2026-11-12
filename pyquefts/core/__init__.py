"""Scientific core of PyQUEFTS (no I/O, no geospatial dependencies)."""

from .conversions import MEHLICH3_TO_OLSEN, P_METHODS, mehlich3_to_olsen, p_to_olsen
from .model import pair_yield, requirements, run_quefts, season_factor, uptake
from .parameters import UPTAKE_ADJUST_DEFAULT, Biomass, CropParameters, Fertilizer, SoilParameters
from .result import R_OUTPUTS, QueftsResult
from .supply import (NutrientSupply, PTOTAL_RATIOS, nut_supply, nut_supply1, nut_supply2,
                     ph_factors_janssen, ph_factors_sattari, temperature_factor_n)

__all__ = [
    "run_quefts", "uptake", "pair_yield", "requirements", "season_factor",
    "SoilParameters", "CropParameters", "Biomass", "Fertilizer", "UPTAKE_ADJUST_DEFAULT",
    "QueftsResult", "R_OUTPUTS",
    "NutrientSupply", "nut_supply", "nut_supply1", "nut_supply2", "PTOTAL_RATIOS",
    "ph_factors_janssen", "ph_factors_sattari", "temperature_factor_n",
    "mehlich3_to_olsen", "p_to_olsen", "P_METHODS", "MEHLICH3_TO_OLSEN",
]
