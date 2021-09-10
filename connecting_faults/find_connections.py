from itertools import product

import geopandas as gpd
from matplotlib import pyplot as plt
import numpy as np
import ezdxf

from fault_mesh.faults import LeapfrogMultiFault
from fault_mesh.connections import ConnectedFaultSystem
from fault_mesh.utilities.graph import connected_nodes
from fault_mesh.utilities.cubit import make_journal_file_multi

data = LeapfrogMultiFault.from_shp("../gis/cfm_gt_1_5.gpkg")

# data = LeapfrogMultiFault.from_shp("../gis/central_cfm.gpkg")
# data = LeapfrogMultiFault.from_shp("../gis/NZ_CFM_v0_9_110821.shp")

dist_tolerance = 200.

connections = []
neighbour_connections = []
for fault in data.faults:
    for other_fault in data.faults:
        if other_fault.name != fault.name:
            if fault.nztm_trace.distance(other_fault.nztm_trace) <= dist_tolerance:
                print(f"Connection: {fault.name} and {other_fault.name}")
                connections.append([fault.name, other_fault.name])
                conditions = [p1.distance(p2) <= dist_tolerance for p1, p2 in product([fault.end1, fault.end2],
                                                                                     [other_fault.end1,
                                                                                      other_fault.end2])]
                if any(conditions):
                    neighbour_connections.append([fault.name, other_fault.name])

major_faults = connected_nodes(neighbour_connections)
data.read_fault_systems("test_central_1_5_composite.csv")
data.generate_curated_faults()
data.read_cutting_hierarchy("test_central_1_5_hierarchy.csv")



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

edges = gpd.GeoDataFrame({"fault_name": [fault.name for fault in data.curated_faults]},
                         geometry=[fault.footprint(smoothed=True) for fault in data.curated_faults])
edges.to_file("edges.gpkg", driver="GPKG")

make_journal_file_multi("edges.gpkg", "footprint_volumes")
edges.plot(ax=ax, alpha=0.5)
ax.set_aspect("equal")
plt.show()
