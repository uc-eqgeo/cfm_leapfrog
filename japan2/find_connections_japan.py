import os.path

import geopandas as gpd
from matplotlib import pyplot as plt
import numpy as np

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.connected import ConnectedFaultSystem

# Read in data
data = LeapfrogMultiFault.from_nz_cfm_shp("example1_faults_trimmed.geojson", depth_type=None, exclude_aus=False)

# Set minimum distance for fault traces to be treated as connected
data.segment_distance_tolerance = 20000.

# Find connections between fault traces
data.find_connections()
data.suggest_fault_systems(out_prefix="example1_combined")

# Read in fault system connections (new file created by editing suggested CSV file)
data.read_fault_systems("example1_combined_suggested.csv")

# Created curated fault system objects (some single LeapfrogFault objects, some ConnectedFaultSystem objects)
data.generate_curated_faults()
data.suggest_cutting_hierarchy(prefix="example1")
# Read cutting hierarchy (in this case unedited)
data.read_cutting_hierarchy("example1_suggested.csv")

# Create directories for output
for dir_name in ["shps", "traces", "end_lines", "footprints"]:
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

# Loop through curated faults and save to shapefiles
plt.close("all")
fig, ax = plt.subplots()
for fault in data.curated_faults:
    # Generate and write depth contours
    fault.generate_depth_contours(np.arange(2000, 32000., 2000.))
    fault.contours.to_file(f"shps/{fault.name}_contours.shp")
    # Write trace
    trace = gpd.GeoSeries(fault.smoothed_trace)
    # Plot trace and depth contours
    trace.plot(ax=ax)
    fault.contours.plot(ax=ax)
    # Write trace to shapefile
    gpd.GeoSeries(fault.smoothed_trace).to_file(f"traces/{fault.name}_trace.shp")
    # Write end lines, although not currently used by leapfrog
    if isinstance(fault, ConnectedFaultSystem):
        gpd.GeoSeries(fault.end_lines(smoothed=True)).to_file(f"end_lines/{fault.name}_end_lines.shp")

# Write footprints to shapefile
footprints = [fault.footprint for fault in data.curated_faults]
for fault in reversed(data.curated_faults):
    fault.adjust_footprint()
    gpd.GeoSeries(fault.footprint).to_file(f"footprints/{fault.name}_footprint.shp")

# Plot footprints
edges = gpd.GeoDataFrame({"fault_name": [fault.name for fault in data.curated_faults]},
                         geometry=[fault.footprint for fault in data.curated_faults])
edges.plot(ax=ax, alpha=0.5)
ax.set_aspect("equal")
plt.show()
