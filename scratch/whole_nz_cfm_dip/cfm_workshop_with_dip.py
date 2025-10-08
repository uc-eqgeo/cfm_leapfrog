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
cfm_shp = r"NZ_CFM_v1_0_rs1km_modified_GDM.gpkg"
# cfm_shp = r"NZ_CFM_v1_0_no_tvz_nztm.gpkg"

# faults = LeapfrogMultiFault.gdf_from_nz_cfm_shp(cfm_shp, exclude_region_polygons=poly_ls)
# faults.crs = 2193
# faults.to_file("NZ_CFM_v1_0_no_tvz_nztm.gpkg", driver="GPKG")


all_faults = LeapfrogMultiFault.from_nz_cfm_shp(cfm_shp, exclude_region_polygons=poly_ls, remove_colons=True)
all_faults.segment_distance_tolerance = 2000.
all_faults.find_connections(verbose=False)
all_faults.read_fault_systems("NZ_CFM_v1_0_rs1km_modified_GDM_connected_edited.csv")
all_faults.generate_curated_faults()

whv = all_faults.curated_fault_dict["Wellington Hutt Valley combined"]
whv4 = whv.segments[3]