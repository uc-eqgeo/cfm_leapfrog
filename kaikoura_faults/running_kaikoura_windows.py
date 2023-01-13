import os.path

import geopandas as gpd
from matplotlib import pyplot as plt
import numpy as np

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.connected import ConnectedFaultSystem
from fault_mesh.utilities.graph import connected_nodes

#read in faults
# note D90 is the seismogenic thickness, whereas Dfc is the theoretical maximum rupture depth (see p. 30 of the CFM report)
# warnings about expected fields because of changes in naming conventions of CFM - should be fine
#easiest way to exclude faults is to edit this gpkg
script_dir=os.path.abspath(".")
data = LeapfrogMultiFault.from_nz_cfm_shp(os.path.join(script_dir,"faults_in_windows","qgis_files","cfm_cut_nonzero.gpkg"), remove_colons=True,
                                        exclude_zero=False, depth_type="Dfc", exclude_aus=False)

# suggest faults to combine (these need editing in csv)
data.find_connections()
major_faults = connected_nodes(data.neighbour_connections)
data.suggest_fault_systems(out_prefix=os.path.join(script_dir,"faults_in_windows","kaikoura"))
#read in edited faults
data.read_fault_systems(fault_system_csv=os.path.join(script_dir,"faults_in_windows","kaikoura_suggested_faults_edited_cfm_v3.csv"),tol=300.)


data.generate_curated_faults()


data.suggest_cutting_hierarchy(os.path.join(script_dir,"faults_in_windows","kaikoura"))
#removing particular fault from hierarchy file
# with open(os.path.join("./faults_in_windows","kaikoura_hierarchy_0slip2rm.csv"),'r') as fin:
#     rm_list = fin.readlines()
# #faults2rm=[fault.strip() for fault in rm_list]
# with open(os.path.join("./faults_in_windows","kaikoura_hierarchy151222.csv"), "r") as hfile:
#     names = hfile.readlines()
# #all_faults = [name.strip() for name in list(names)]
# new_hierarchy=[fault for fault in names if not fault in rm_list]
# with open(os.path.join("./faults_in_windows","kaikoura_hierarchy151222_nonzero.csv"), "w") as fout:
#     fout.writelines(new_hierarchy)


data.read_cutting_hierarchy(os.path.join(script_dir,"faults_in_windows","kaikoura_hierarchy151222_nonzero.csv"))

for dir_name in ["shps", "traces", "end_lines", "footprints"]:
    if not os.path.exists(os.path.join(script_dir,"faults_in_windows",dir_name)):
        os.mkdir(os.path.join(script_dir,"faults_in_windows",dir_name))




plt.close("all")
fig, ax = plt.subplots()
for fault in data.curated_faults:
    print(fault.name)


    fault.generate_depth_contours(np.arange(2000, 32000., 2000.))


    fault.contours.to_file(os.path.join(script_dir,"faults_in_windows",f"shps/{fault.name}_contours.shp"))
    trace = gpd.GeoSeries(fault.smoothed_trace)
    trace.plot(ax=ax)
    fault.contours.plot(ax=ax)

    gpd.GeoSeries(fault.smoothed_trace).to_file(os.path.join(script_dir,"faults_in_windows",f"traces/{fault.name}_trace.shp"))
    if isinstance(fault, ConnectedFaultSystem):
        print("connected")
        gpd.GeoSeries(fault.end_lines(smoothed=True)).to_file(os.path.join(script_dir,"faults_in_windows",f"end_lines/{fault.name}_end_lines.shp"))

footprints = [fault.footprint for fault in data.curated_faults]
for fault in reversed(data.curated_faults):
    fault.adjust_footprint()
    gpd.GeoSeries(fault.footprint).to_file(os.path.join(script_dir,"faults_in_windows",f"footprints/{fault.name}_footprint.shp"))

edges = gpd.GeoDataFrame({"fault_name": [fault.name for fault in data.curated_faults]},
                         geometry=[fault.footprint for fault in data.curated_faults])
edges.to_file(os.path.join(script_dir,"faults_in_windows","edges.gpkg"), driver="GPKG")
#
# make_journal_file_multi("edges.gpkg", "footprint_volumes")
edges.plot(ax=ax, alpha=0.5)
ax.set_aspect("equal")
plt.show()
