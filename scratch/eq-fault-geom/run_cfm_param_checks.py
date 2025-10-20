import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

from fault_mesh.faults.opensha_faults import CfmMultiFault
from fault_mesh.io.shapefile_operations import geopandas_polygon_to_gmt

shp = "data/cfm1_0.gpkg"
# shp = "../../../data/cfm_shapefile/cfm_0_9_stirling_depths.gpkg"

# Polygons to exclude faults from XML
exclude_shp = "data/bop_exclusion.gpkg"
exclude_df = gpd.GeoDataFrame.from_file(exclude_shp)
# Change to WGS for use with Matt's TVZ polygon
exclude_df_wgs = exclude_df.to_crs(epsg=4326)
poly_ls = list(exclude_df_wgs.geometry.explode())

# To read in Matt's TVZ polygon
matt_array = np.loadtxt("data/tvz_polygon.txt")

# Polygon requires lon lat (rather than lat lon)
matt_poly = Polygon(matt_array)
poly_ls.append(matt_poly)

# # Write to shapefile for visualization in GIS
matt_gs = gpd.GeoSeries(matt_poly, crs=4326)
matt_gs.to_file("matt_poly.shp")

# Width of trace buffer (in metres)
buffer_width = 5000.

# Read and write data
data_all = CfmMultiFault.from_shp(shp, depth_type=f"Dfc", min_dfc=3.9)
data_notvz = CfmMultiFault.from_shp(shp, exclude_region_polygons=poly_ls,
                                        exclude_region_min_sr=1.8, depth_type=f"Dfc", min_dfc=3.9)


# Fault polygons for any non-horizontal faults
# polygons = [fault.combined_buffer_polygon(buffer_width) for fault in data_d90_notvz.faults if abs(fault.down_dip_vector[-1]) > 1.e-3]
# Write out polygons
# polygons_gdf = gpd.GeoDataFrame(geometry=polygons, crs=4326)
# polygons_gdf.to_file("fault_buffers.shp")
# Write out for plotting in matlab
# geopandas_polygon_to_gmt(polygons_gdf, "fault_polygons_notvz.txt", matlab=True)

# Write out XML files
for file_handle, dataset in zip(["all", "no_tvz"],
                               [data_all, data_notvz]):
    xml_buffer = dataset.to_opensha_xml(exclude_subduction=True, buffer_width=buffer_width, write_buffers=False)
    with open(f"cfm_1_0A_{file_handle}.xml", "wb") as f:
        f.write(xml_buffer)
