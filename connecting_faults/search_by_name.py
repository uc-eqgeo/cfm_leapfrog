import geopandas as gpd
import fnmatch
from fault_mesh.faults import LeapfrogMultiFault

data_d90_all = LeapfrogMultiFault.from_shp("../gis/cfm_gt_1_5.gpkg", depth_type="D90")


for fault in data_d90_all.faults:
    if fnmatch.fnmatch(fault.name, "Alpine*"):
        print(fault.name)


