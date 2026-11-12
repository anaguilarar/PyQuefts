"""Module 1: run QUEFTS from user-supplied soil values (mirrors the reagro.org tutorial)."""

import numpy as np

import pyquefts as pq
from pyquefts import analysis as an

# --- Step 1: native soil supply from topsoil properties (0-20 cm)
sup = pq.nut_supply2(temp=25, ph=6.0, soc=15, kex=1.0, polsen=10, ptotal=500)
print("Base supply (kg/ha):", {k: round(float(v), 2) for k, v in sup.to_dict().items()})

# --- Steps 2-4: a barley crop, 100 kg N/ha, 110-day season
q = pq.Quefts(soil=sup.to_soil(N_recovery=0.5, P_recovery=0.1, K_recovery=0.5), crop="barley",
              fert=dict(N=100, P=0, K=0),
              biom=dict(leaf_att=2200, stem_att=2700, store_att=4800, SeasonLength=110))
res = q.run()
print(res)
print("Intermediate yield ranges (kg/ha):", {k: round(res.scalar(k)) for k in ("YNA", "YND", "YPA", "YPD", "YKA", "YKD")})
print("Most limiting nutrient:", an.limiting_nutrient(res))

# --- one-call version from raw soil values
res2 = pq.run_point(ph=6.0, soc=15, kex=1.0, polsen=10, ptotal=500, temp=25, crop="barley",
                    fert=dict(N=100, P=0, K=0), biom=dict(leaf_att=2200, stem_att=2700, store_att=4800, SeasonLength=110))
assert res2.scalar("store_lim") == res.scalar("store_lim")

# --- N response curve (vectorized: one model call for all rates)
curve = an.response_curve("N", np.arange(0, 201, 25), "rice", fert=dict(N=0, P=0, K=0),
                          biom=dict(leaf_att=2651, stem_att=5053, store_att=8208))
print(curve.round(1).to_string(index=False))

# --- N x P response surface at K = 0 and 200
surf = an.response_surface(np.arange(0, 401, 50), np.arange(0, 401, 50), [0, 200], "rice",
                           biom=dict(leaf_att=2651, stem_att=5053, store_att=8208))
print(surf.sel(K=200).round(0).to_pandas())

# --- balanced nutrient uptake requirements (Yang et al. 2017 style)
print(an.balanced_uptake("maize", [2000, 4000, 6000, 8000], ymax=10000, leaf_ratio=0.46, stem_ratio=0.56).round(1))
