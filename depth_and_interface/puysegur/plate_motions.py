import geopandas as gpd
import numpy as np
from euler_pole import EulerPole
from pyproj import Transformer
from shapely.geometry import Point

wgs2nztm = Transformer.from_crs("epsg:4326", "epsg:2193", always_xy=True)
nztm2wgs = Transformer.from_crs("epsg:2193", "epsg:4326", always_xy=True)

border = gpd.read_file("nz_region.geojson")

x1, y1, x2, y2 = border.total_bounds

x = np.arange(x1, x2, 100000.)
y = np.arange(y1, y2, 100000.)

xv, yv = np.meshgrid(x, y)

x_flat = xv.flatten()
y_flat = yv.flatten()

lon, lat = nztm2wgs.transform(x_flat, y_flat)
aust_paci = EulerPole(-61.04, 184.19 - 360., 1.08)

azimuths = []
velocities = []
for lon_i, lat_i in zip(lon, lat):
    az_i, vel_i = aust_paci.velocity(lat_i, lon_i)
    azimuths.append(az_i)
    velocities.append(vel_i)

velocities = np.array(velocities)
azimuths = np.array(azimuths)

velocity_gdf = gpd.GeoDataFrame({"velocity": velocities, "azimuth": azimuths},
                                geometry=[Point(x, y) for x, y in zip(x_flat, y_flat)], crs=2193)
velocity_gdf.to_file("nz_vels.geojson", driver="GeoJSON")


