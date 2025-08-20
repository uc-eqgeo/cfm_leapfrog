from fault_mesh.utilities.meshing import triangulate_contours
from fault_mesh.utilities.splines import spline_fit_contours
import geopandas as gpd
from matplotlib import pyplot as plt

# Import modules 
from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import os
import numpy as np
import geopandas as gpd

# Set coordinate system (optional) EPSG code
# If not desired, set to None
epsg = 2193

# Read in fault data from shapefile
fault_fname = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\NZ_CFM_v1_0_rs1km_modified_GDM.gpkg"
fault_data = LeapfrogMultiFault.from_nz_cfm_shp(fault_fname, remove_colons=True, epsg=epsg, smoothing_n=None, dip_multiplier=1.0)

dist_tolerance = 200.
fault_data.find_connections(verbose=False)

fault_data.read_fault_systems(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\NZ_CFM_v1_0_rs1km_modified_GDM_connected_edited.csv")
fault_data.generate_curated_faults()

alpine = fault_data.name_dict['Alpine George - Wairau combined']
alpine.generate_depth_contours(np.arange(0., 32000., 1000.), smoothing=False)

spline_contours = spline_fit_contours(alpine.contours, point_spacing=100., output_spacing=1000.)

spline_contours.plot()
plt.show()

for fault in fault_data.connected_faults:
    print(fault.name)
    fault.generate_depth_contours(np.arange(0., 32000., 1000.), smoothing=False)
    spline_contours = spline_fit_contours(fault.contours, point_spacing=100., output_spacing=1000.)
    contour_file = os.path.join(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_contours", f"{fault.name}_depth_contours.geojson")
    spline_contours.to_file(contour_file, driver="GeoJSON")
    vtk_file = os.path.join(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks", f"{fault.name}_depth_contours.vtk")
    triangulate_contours(spline_contours, vtk_file, mesh_format="vtk", check_mesh=False)