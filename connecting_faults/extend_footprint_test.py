import numpy as np

from fault_mesh.faults import LeapfrogMultiFault
from fault_mesh.utilities.graph import connected_nodes

data = LeapfrogMultiFault.from_shp("../gis/cfm_gt_1_5.gpkg")

# data = LeapfrogMultiFault.from_shp("../gis/central_cfm.gpkg")
# data = LeapfrogMultiFault.from_shp("../gis/NZ_CFM_v0_9_110821.shp")

dist_tolerance = 200.

data.find_connections()

major_faults = connected_nodes(data.neighbour_connections)
data.read_fault_systems("test_central_1_5_composite.csv")
data.generate_curated_faults()
data.read_cutting_hierarchy("test_central_1_5_hierarchy.csv")

pap = data.curated_fault_dict["Papatea"]

jkn = data.curated_fault_dict["Jordan - Kek - Needles combined"]
hope = data.curated_fault_dict["Hope combined"]

closest_jkn = min(jkn.segments, key=lambda x: x.nztm_trace.distance(hope.nztm_trace))
closest_hope = min(hope.segments, key=lambda x: x.nztm_trace.distance(jkn.nztm_trace))

closest_jkn.generate_depth_contours(np.arange(2000, 32000., 2000.))
hope.generate_depth_contours(np.arange(2000, 32000., 2000.))
