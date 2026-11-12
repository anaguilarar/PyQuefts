"""Module 2: simulate at coordinates using the Honduras soil data cube."""

from pyquefts.geo import SoilCube, extract_points, run_points, sample_at_points, season_mean_temperature

CUBE = r"D:/OneDrive - CGIAR/projects/spatial_files/raster/soil/hnd_sh_soilgrids_format.nc"

# the Suelos de Honduras P is Mehlich-3: convert to Olsen P (McDowell et al. 2023, non-calcareous)
cube = SoilCube.open(CUBE, chunks={"depth": -1, "y": 1024, "x": 1024}, p_method="mehlich3")
print(cube.describe())

lon = [-87.65, -86.20, -88.90]
lat = [14.45, 14.80, 15.50]

# soil inputs in model units for the 0-20 cm layer (thickness-weighted over 0-5, 5-15, 15-30 cm)
print(extract_points(cube, lon, lat).round(2).to_string())

# growing-season mean temperature from the CHIRTS-ERA5 daily cube (primera season, 1991-2020 climatology)
temp = season_mean_temperature(start="05-15", end="08-31", years=range(1991, 2021), bbox=(-89.4, 12.9, -83.1, 16.6))
print("season temperature at points:", sample_at_points(temp, lon, lat).round(1))

# full run: Sattari (2014) supply, maize, 6 t/ha attainable grain, no fertilizer
df = run_points(cube, lon, lat, crop="maize", temp=temp, store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
                fert=dict(N=0, P=0, K=0))
cols = ["lon", "lat", "temp", "ph", "soc", "kex", "polsen", "N_base_supply", "P_base_supply", "K_base_supply",
        "store_lim", "N_uptake", "P_uptake", "K_uptake", "N_gap", "P_gap", "K_gap"]
print(df[cols].round(1).to_string())

# --- N response curve at each point (soil extracted once; all points x all rates in one vectorized run)
import numpy as np  # noqa: E402
from pyquefts.geo import response_points  # noqa: E402

resp = response_points(cube, lon, lat, crop="maize", nutrient="N", rates=np.arange(0, 201, 20), temp=temp,
                       store_att=6000, leaf_ratio=0.46, stem_ratio=0.56, fert=dict(N=0, P=0, K=0))
print(resp.pivot(index="N", columns="label", values="store_lim").round(0).to_string())
try:
    import matplotlib
    matplotlib.use("Agg")
    from pyquefts import plot_curves
    axes = plot_curves(resp, y=["store_lim", "N_gap"])
    axes[0].figure.tight_layout()
    axes[0].figure.savefig("points_N_response.png", dpi=120)
    print("wrote points_N_response.png")
except ImportError:
    print("matplotlib not installed: skipping the plot")
