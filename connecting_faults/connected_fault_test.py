from faults.connected import ConnectedFaultSystem
from fault_mesh.faults import LeapfrogMultiFault
import geopandas as gpd

# data_d90_all = LeapfrogMultiFault.from_shp("../gis/cfm_gt_1_5.gpkg")
data_d90_all = LeapfrogMultiFault.from_shp("../gis/central_cfm.gpkg")

test = ConnectedFaultSystem("test", ["Hope*", "Kelly"], data_d90_all, excluded_names=["*Te Rapa*", "*Hanmer NW", "*Tarama*"],
                            smooth_trace_refinements=5)

smoothed_trace = gpd.GeoSeries(test.smoothed_overall_trace)
proj_points = gpd.GeoSeries(test.segment_boundaries)
smoothed_segments = gpd.GeoSeries(test.smoothed_segments[:2])

# plt.close("all")
# fig, ax = plt.subplots()
# smoothed_trace.plot(ax=ax)
# proj_points.plot(ax=ax)
# smoothed_segments.plot(ax=ax, color="r")
#
# plt.show()


