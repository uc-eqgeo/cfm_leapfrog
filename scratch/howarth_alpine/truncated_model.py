# Import modules 
from fault_mesh.faults.leapfrog import LeapfrogMultiFault, ConnectedFaultSystem
import os
import numpy as np
import geopandas as gpd

# Set coordinate system (optional) EPSG code
# If not desired, set to None
epsg = 2193
# Read in fault data from shapefile
fault_data = LeapfrogMultiFault.from_nz_cfm_shp("alpine_cfm_emily_truncated_1km.geojson", remove_colons=True, epsg=epsg, smoothing_n=None)
dist_tolerance = 200.0  # meters

fault_data.find_connections(verbose=False)
fault_data.suggest_fault_systems("alpine_connected")

fault_data.read_fault_systems("alpine_connected_truncated_edited.csv")
fault_data.generate_curated_faults()

fault_data.suggest_cutting_hierarchy("alpine_hierarchy")
fault_data.read_cutting_hierarchy("alpine_hierarchy_suggested.csv")

# for dir_name in ["depth_contours", "traces", "footprints", "footprints_lines"]:
#     if not os.path.exists(dir_name):
#         os.mkdir(dir_name)

for fault in fault_data.curated_faults:
    # Generate depth contours
    depths = np.arange(0., 20000., 500.)
    fault.generate_depth_contours(depths, smoothing=False)
    # Write contours to files
    fault.contours.to_file(f"depth_contours/{fault.name}_truncated_contours.shp")
    # Write traces
    fault.nztm_trace_geoseries.to_file(f"traces/{fault.name}_truncated_trace.shp")
    # write endpoints
    if isinstance(fault, ConnectedFaultSystem):
        end_points = fault.generate_sr_rake_points(depths, smoothing=False)
        end_points.to_file(f"endpoints/{fault.name}_truncated_dip_endpoints.shp")

# Write fault footprints
for fault in reversed(fault_data.curated_faults):
    fault.adjust_footprint()
    fault.footprint_geoseries.to_file(f"footprints/{fault.name}_footprint.shp")
    # Newer versions of Leapfrog Geo require fault footprints to be written out as lines
    fault.footprint_geoseries.boundary.to_file(f"footprints_lines/{fault.name}_truncated_footprint.shp")