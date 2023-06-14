import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString, MultiLineString, MultiPolygon, Point, MultiPoint
from shapely.ops import unary_union
import meshio
import numpy as np
from glob import glob
import os.path

# Extract data of interest from Ellis et al. (2022) and save to CSV.
ellis_data = pd.read_csv("ESupp1_RuptureDepthSurfaces_v1.0.csv")
ellis_data["Dfcomb(m)"] = ellis_data["Dfcomb(km)"] * -1000
trimmed_data = ellis_data[["Easting(NZTM)", "Northing(NZTM)", "Dfcomb(m)"]]
trimmed_data.to_csv("ellis_dfcomb_surface.csv", index=False)

# Get bounds of CFM model
cfm = gpd.read_file("NZ_CFM_v1_0_shapefile/NZ_CFM_v1_0.shp")
cfm_bounds = cfm.total_bounds
# Add buffer to bounds
buffer = 500000.
cfm_bounds_gmt = cfm_bounds + np.array([-buffer, -buffer, buffer, buffer])
cfm_bounds_gmt = [cfm_bounds_gmt[0], cfm_bounds_gmt[2], cfm_bounds_gmt[1], cfm_bounds_gmt[3]]
# Print bounds in GMT format (not actually used, would need to add rounding to be used directly at present)
print("bounds=-R{}/{}/{}/{}".format(*cfm_bounds_gmt))

buffer_distance = 2.5e4
# Read in mesh, get triangle vertices and outline to use in gridding
mesh = meshio.read("hik_kerm3k.stl")
mesh_points = mesh.points
mesh_cells = mesh.cells_dict["triangle"]
mesh_tri = mesh_points[mesh_cells]
mesh_poly = gpd.GeoDataFrame(geometry=gpd.GeoSeries([Polygon(tri) for tri in mesh_tri], crs=2193))
# This mesh didn't have the southern cut off, so we need to trim it
# outline = unary_union(mesh_poly.geometry.values)
# outline_gs = gpd.GeoSeries(outline, crs=2193)
# outline_gs.to_file("hikurangi_outline.geojson", driver="GeoJSON")
# Read in trimmed outline after editing in QGIS
outline_gs = gpd.read_file("hikurangi_outline_trimmed.geojson")
outline = list(outline_gs.geometry.values[0].geoms)[0]
#
outline_array = np.array(outline.exterior.coords)
buffered_hikurangi = outline.buffer(buffer_distance)
np.savetxt("hikurangi_points.csv", mesh_points, delimiter=" ")
np.savetxt("hikurangi_outline.csv", outline_array, delimiter=" ")
np.savetxt("hikurangi_outline_buffered.csv", np.array(buffered_hikurangi.exterior.coords), delimiter=" ")

puysegur_points = np.loadtxt("puysegur_points.txt")
puysegur_outline = gpd.read_file("puysegur_outline.geojson")
puysegur_buffered = puysegur_outline.geometry.values[0].buffer(buffer_distance)
puysegur_outline_array = np.array(puysegur_outline.geometry.values[0].exterior.coords)
np.savetxt("puysegur_outline.gmt", puysegur_outline_array, delimiter=" ")
np.savetxt("puysegur_outline_buffered.gmt", np.array(puysegur_buffered.exterior.coords), delimiter=" ")

# Hannu's alterations
extra_bits = list(glob("Leapfrog/*.shp"))
for bit in extra_bits:
    fname = os.path.basename(bit).split(".")[0]
    coords = np.array(gpd.read_file(bit).geometry.values[0].exterior.coords)
    np.savetxt(f"{fname}.gmt", coords, delimiter=" ")





# mesh_tri = mesh_tri.reshape((mesh_tri.shape[0], 9))S