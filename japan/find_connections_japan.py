import os.path

import geopandas as gpd
from matplotlib import pyplot as plt
import numpy as np

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.connected import ConnectedFaultSystem
from fault_mesh.utilities.graph import connected_nodes

data = LeapfrogMultiFault.from_nz_cfm_shp("example1_faults_trimmed.geojson", depth_type=None)


data.segment_distance_tolerance = 1000.

# data.find_connections()
# data.suggest_fault_systems("example1")

data.read_fault_systems("example1_combined_edited.csv")
data.generate_curated_faults()
data.suggest_cutting_hierarchy("example1")
data.read_cutting_hierarchy("example1_suggested.csv")
#
for dir_name in ["shps", "traces", "end_lines", "footprints"]:
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
#
plt.close("all")
fig, ax = plt.subplots()
for fault in data.curated_faults:
    fault.generate_depth_contours(np.arange(2000, 32000., 2000.))
    fault.contours.to_file(f"shps/{fault.name}_contours.shp")
    trace = gpd.GeoSeries(fault.smoothed_trace)
    trace.plot(ax=ax)
    fault.contours.plot(ax=ax)

    gpd.GeoSeries(fault.smoothed_trace).to_file(f"traces/{fault.name}_trace.shp")
    if isinstance(fault, ConnectedFaultSystem):
        gpd.GeoSeries(fault.end_lines(smoothed=True)).to_file(f"end_lines/{fault.name}_end_lines.shp")

footprints = [fault.footprint for fault in data.curated_faults]
for fault in reversed(data.curated_faults):
    fault.adjust_footprint()
    gpd.GeoSeries(fault.footprint).to_file(f"footprints/{fault.name}_footprint.shp")

edges = gpd.GeoDataFrame({"fault_name": [fault.name for fault in data.curated_faults]},
                         geometry=[fault.footprint for fault in data.curated_faults])
edges.to_file("edges.gpkg", driver="GPKG")
#
# make_journal_file_multi("edges.gpkg", "footprint_volumes")
edges.plot(ax=ax, alpha=0.5)
ax.set_aspect("equal")
plt.show()