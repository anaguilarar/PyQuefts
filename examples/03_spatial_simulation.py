"""Module 3: spatial QUEFTS run over a region of the Honduras soil cube, written to NetCDF and GeoTIFF."""

import time

from pyquefts.geo import SoilCube, run_spatial, season_mean_temperature, write_geotiff, write_netcdf

CUBE = r"D:/OneDrive - CGIAR/projects/spatial_files/raster/soil/hnd_sh_soilgrids_format.nc"
BBOX = (-88.0, 14.0, -87.5, 14.5)  # xmin, ymin, xmax, ymax (degrees)

# the Suelos de Honduras P is Mehlich-3: convert to Olsen P (McDowell et al. 2023, non-calcareous)
cube = SoilCube.open(CUBE, chunks={"depth": -1, "y": 1024, "x": 1024}, p_method="mehlich3")

# growing-season temperature (primera, 1991-2020) from CHIRTS-ERA5, regridded onto the soil grid inside run_spatial
temp = season_mean_temperature(start="05-15", end="08-31", years=range(1991, 2021), bbox=BBOX)

t = time.time()
ds = run_spatial(
    cube, "maize", store_att=6000, leaf_ratio=0.46, stem_ratio=0.56,
    fert=dict(N=60, P=15, K=30), temp=temp, supply_method="sattari2014", ptotal_ratio="afsis",
    bbox=BBOX, top=0, bottom=20, compute=True,
)
print(ds)
print(f"{ds.store_lim.size:,} cells in {time.time() - t:.1f} s")
print("nutrient-limited grain yield (kg/ha DM): ",
      float(ds.store_lim.min()), float(ds.store_lim.mean()), float(ds.store_lim.max()))

write_netcdf(ds, "quefts_maize_comayagua.nc")
write_geotiff(ds, "quefts_maize_comayagua.tif", variables=["store_lim", "N_gap", "P_gap", "K_gap"])

# A scenario comparison stays lazy until computed; e.g. yield gain from 60 kg N/ha
ds0 = run_spatial(cube, "maize", store_att=6000, leaf_ratio=0.46, stem_ratio=0.56, fert=dict(N=0, P=15, K=30),
                  temp=temp, bbox=BBOX, keep_inputs=False)
gain = (ds.store_lim - ds0.store_lim).compute()
print("mean yield gain from 60 kg N/ha:", float(gain.mean()))
