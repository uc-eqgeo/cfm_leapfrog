from itertools import product

import geopandas as gpd
from matplotlib import pyplot as plt

from fault_mesh.faults import LeapfrogMultiFault
from fault_mesh.utilities.graph import connected_nodes

data = LeapfrogMultiFault.from_shp("../gis/cfm_gt_1_5.gpkg")

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
