# Import modules 
from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import os
import numpy as np
import geopandas as gpd

# Set coordinate system (optional) EPSG code
# If not desired, set to None
epsg = 2193
# Read in fault data from shapefile
fault_data = LeapfrogMultiFault.from_nz_cfm_shp("alpine_cfm_emily_multisegment_1km.geojson", remove_colons=True, epsg=epsg, smoothing_n=None)
dist_tolerance = 200.0  # meters

for fault in fault_data.curated_faults:
    # Generate depth contours
    fault.generate_depth_contours(np.arange(0., 22000., 2000.), smoothing=False)
    fault.contours.to_file(f"depth_contours/multisegment/{fault.name}_planar_contours.shp")
