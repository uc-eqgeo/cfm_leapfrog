import geopandas as gpd
from fault_mesh.utilities.merging import merge_multiple_nearly_adjacent_segments
from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import matplotlib.pyplot as plt
from shapely.geometry import LineString, MultiPoint
import numpy as np

data = LeapfrogMultiFault.from_nz_cfm_shp("NZ_CFM_v1_0_rs1km_GDM_1mmyr_modified_Alpine_Caswell.gpkg", remove_colons=True,
                                            exclude_aus=False, exclude_zero=False, exclude_region_min_sr=0.0,
                                          smoothing_n=None, epsg=2193)

data.read_fault_systems("NZ_CFM_v1_0_rs1km_GDM_1mmyr_modified_Alpine_Caswell_connected_edited.csv")
f0 = data.connected_faults[0]
f0.generate_depth_contours(np.arange(2000, 32000., 2000.))
f0.contours.to_file("contours.shp")

f0.nztm_trace_geoseries.to_file("trace.shp")
