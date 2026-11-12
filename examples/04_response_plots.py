"""Plotting extension: nutrient response curves, surfaces and the uptake-yield envelope.

Requires matplotlib (``pip install pyquefts[plot]``). Writes PNG files to the
current directory.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import pyquefts as pq  # noqa: E402
from pyquefts import analysis as an  # noqa: E402

# soil supply from topsoil properties (Olsen P here; pass p_method="mehlich3" in run_point / SoilCube for Mehlich-3)
sup = pq.nut_supply2(temp=26, ph=6.0, soc=15, kex=1.0, polsen=10)

# --- 1. N response: the Rquefts-tutorial loop, vectorized. Varying one input in `fert` gives array outputs
#        and result.plot() puts that input on the x axis.
q = pq.Quefts(soil=sup.to_soil(), crop="maize",
              fert=dict(N=np.arange(0, 160, 10), P=10, K=60),
              biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
res = q.run()
print(res.to_dataframe(include_inputs=True)[["N", "store_lim", "N_uptake", "N_gap"]].round(1).to_string(index=False))
axes = res.plot(y=["store_lim", "N_uptake", "N_gap"])
axes[0].figure.suptitle("Maize: response to N (P = 10, K = 60 kg/ha)")
axes[0].figure.tight_layout()
axes[0].figure.savefig("response_N.png", dpi=130)

# --- 2. N, P and K panels in one call (other two nutrients held at q.fert)
q.set(N=60)
axes = q.plot(rates=np.arange(0, 201, 10))
axes[0].figure.tight_layout()
axes[0].figure.savefig("response_NPK.png", dpi=130)

# --- 3. N response at several attainable yields (water-limited yield levels)
N = np.arange(0, 201, 10)[:, None]
store = np.array([3000.0, 5000.0, 8000.0])[None, :]
q2 = pq.Quefts(soil=sup.to_soil(), crop="maize", fert=dict(N=N, P=20, K=60),
               biom=dict(leaf_att=0.46 * store, stem_att=0.56 * store, store_att=store, SeasonLength=120))
ax = q2.run().plot(x="N", hue="store_att")
ax.set_title("Maize: N response at three attainable yields")
ax.figure.tight_layout()
ax.figure.savefig("response_N_by_yatt.png", dpi=130)

# --- 4. N x P response surface at K = 60
rs = an.response_surface(np.arange(0, 201, 20), np.arange(0, 61, 5), [60], "maize", soil=sup.to_soil(),
                         biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
ax = pq.plot_surface(rs, x="N", y="P")
ax.figure.tight_layout()
ax.figure.savefig("surface_NP.png", dpi=130)

# --- 5. uptake-yield envelope (Janssen 1990 Fig. 1) with the runs of step 1 overlaid
env = an.uptake_yield_envelope("maize", biom=dict(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120))
axes = pq.plot_envelope(env, res)
axes[0].figure.suptitle("Maize: uptake-yield envelope and modelled runs")
axes[0].figure.tight_layout()
axes[0].figure.savefig("envelope.png", dpi=130)
plt.close("all")
print("wrote response_N.png, response_NPK.png, response_N_by_yatt.png, surface_NP.png, envelope.png")
