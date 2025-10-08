from fault_mesh.utilities.meshing import triangulate_contours
from fault_mesh.utilities.splines import spline_fit_contours
import geopandas as gpd
from matplotlib import pyplot as plt
import meshio

# Import modules 
from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import os
import numpy as np
import geopandas as gpd
from shapely import LineString, MultiLineString

# Set coordinate system (optional) EPSG code
# If not desired, set to None
epsg = 2193

# Read in fault data from shapefile
fault_fname = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\NZ_CFM_v1_0_rs1km_modified_GDM.gpkg"

dist_tolerance = 1000.
fault_data = LeapfrogMultiFault.from_nz_cfm_shp(fault_fname, remove_colons=True, epsg=epsg, smoothing_n=None, dip_multiplier=1.0)
fault_data.segment_distance_tolerance = dist_tolerance
fault_data.find_connections(verbose=False)

fault_data.read_fault_systems(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\NZ_CFM_v1_0_rs1km_modified_GDM_connected_edited.csv")
fault_data.generate_curated_faults()

alpine = fault_data.name_dict['Alpine George - Wairau combined']
alpine.generate_depth_contours(np.arange(0., 32000., 1000.), smoothing=False)

fault_data.read_cutting_hierarchy(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\NZ_CFM_v1_0_rs1km_modified_GDM_Final\NZ_CFM_v1_0_rs1km_modified_GDM_hierarchy_edited.csv")

# spline_contours = spline_fit_contours(alpine.contours, point_spacing=100., output_spacing=1000.)

# hope = fault_data.name_dict['Hope combined']
fault_data.adjust_traces_for_terminations(fit_distance=5.e3, extend_distance=40.e3, resolution=1.e3)

adjusted = gpd.GeoDataFrame({"name": [fault.name for fault in fault_data.curated_faults]}, geometry=[fault.nztm_trace for fault in fault_data.curated_faults], crs=epsg)
adjusted.to_file("test.geojson", driver="GeoJSON")

waimana = fault_data.name_dict['Waimana combined']
waimana.adjust_trace()

leader = fault_data.name_dict['Leader combined']

# spline_contours.plot()
# plt.show()

for fault in fault_data.curated_faults:
    print(fault.name)
    try:
        fault.generate_depth_contours(np.arange(0., 32000., 2000.), smoothing=False)
        mesh = fault.mesh_fault_surface(check_mesh=False)
    except Exception as e:
        print(f"Failed to mesh {fault.name}: {e}")
        fault.contours.to_file(f"failed_contours_{fault.name}.geojson", driver="GeoJSON")
        mesh = fault.mesh_simple_contours(np.arange(0., 32000., 2000.))
    vtk_file = os.path.join(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks", f"{fault.name}_depth_contours.vtk")
    geojson_file = os.path.join(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_contours", f"{fault.name}_depth_contours.geojson")
    fault.contours.to_file(geojson_file, driver="GeoJSON")
    mesh.write(vtk_file)

    # spline_contours = spline_fit_contours(fault.contours, point_spacing=100., output_spacing=1000.)
    # contour_file = os.path.join(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_contours", f"{fault.name}_depth_contours.geojson")
    # spline_contours.to_file(contour_file, driver="GeoJSON")
    
    # triangulate_contours(spline_contours, vtk_file, mesh_format="vtk", check_mesh=False)


contours_list = []
for contour in alpine.contours.geometry.values:
    if isinstance(contour, LineString):
        contours_list.append(np.array(contour.coords))
    else:
        assert isinstance(contour, MultiLineString)
        contours_list.append(np.vstack([np.array(line.coords) for line in contour.geoms]))


all_contour_points = np.vstack(contours_list)

mesh = meshio.read(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\Alpine George - Wairau combined_depth_contours.vtk")