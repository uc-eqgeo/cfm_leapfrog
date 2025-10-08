from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

# Polygons to exclude faults from XML
exclude_shp = "bop_exclusion.gpkg"
exclude_df = gpd.GeoDataFrame.from_file(exclude_shp)
# Change to WGS for use with Matt's TVZ polygon
exclude_df_wgs = exclude_df.to_crs(epsg=4326)
poly_ls = list(exclude_df_wgs.geometry.explode())

# To read in Matt's TVZ polygon
matt_array = np.loadtxt("tvz_polygon.txt")

# Polygon requires lon lat (rather than lat lon)
matt_poly = Polygon(matt_array)
poly_ls.append(matt_poly)
cfm_shp = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\gis\NZ_CFM_v1_0_shapefile\NZ_CFM_v1_0.shp"

faults = LeapfrogMultiFault.gdf_from_nz_cfm_shp(cfm_shp, exclude_region_polygons=poly_ls)
faults.crs = 2193
faults.to_file("NZ_CFM_v1_0_no_tvz_nztm.gpkg", driver="GPKG")