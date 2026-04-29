import os
import numpy as np
from fault_mesh.faults.leapfrog import LeapfrogMultiFault

# Set output coordinate system EPSG code (set to None to disable reprojection)
epsg = 2193

# Read fault data from GeoJSON file
fault_fname = r"C:\Users\arh128\vscode_projects\cfm_leapfrog\docs\tutorials\tutorial_gis\central_nz_minimal_data.shp"

# Build fault network: set segment connection tolerance and find connections between segments
dist_tolerance = 1000.
fault_data = LeapfrogMultiFault.from_nz_cfm_shp(fault_fname, remove_colons=True, epsg=epsg, smoothing_n=None, dip_multiplier=1.0, exclude_zero=False)
fault_data.segment_distance_tolerance = dist_tolerance
fault_data.find_connections(verbose=False)

# Read fault system definitions and generate curated (multi-segment) faults
fault_data.read_fault_systems(r"C:\Users\arh128\vscode_projects\cfm_leapfrog\docs\tutorials\define_connections_data\central_gt1_5_connected_edited.csv")
fault_data.generate_curated_faults()

# Read fault cutting hierarchy to determine termination relationships
fault_data.read_cutting_hierarchy(r"C:\Users\arh128\vscode_projects\cfm_leapfrog\docs\tutorials\define_connections_data\central_gt1_5_hierarchy_edited.csv")

fault_data.suggest_trace_extensions(out_file="suggested_trace_extensions.csv", fit_distance=5.e3, extend_distance=40.e3, proximity_threshold=1.e3)
fault_data.read_trace_extensions("trace_extensions_edited.csv")
fault_data.apply_trace_extensions()
# fault_data.adjust_traces_for_terminations(fit_distance=5.e3, extend_distance=40.e3, resolution=1.e3)

# Create output directories for meshes and contours
if not os.path.exists(r"test_vtks"):
    os.mkdir(r"test_vtks")
if not os.path.exists(r"test_contours"):
    os.mkdir(r"test_contours")
if not os.path.exists(r"test_objs"):
    os.mkdir(r"test_objs")

# Generate depth contours and triangulated mesh for each fault surface
for fault in fault_data.curated_faults:
    print(fault.name)
    try:
        fault.generate_depth_contours(np.arange(0., 32000., 500.), smoothing=False)
        mesh = fault.mesh_fault_surface(check_mesh=False, resolution=500.)
    except Exception as e:
        # On meshing failure, fall back to simple contour meshing and save failed contours for inspection
        print(f"Failed to mesh {fault.name}: {e}")
        fault.contours.to_file(f"failed_contours_{fault.name}.geojson", driver="GeoJSON")
        mesh = fault.mesh_simple_contours(np.arange(0., 32000., 500.))
    #Write depth contours as GeoJSON and mesh as VTK and OBJ for inspection
    vtk_file = os.path.join(r"test_vtks", f"{fault.name}_depth_contours.vtk")
    obj_file = os.path.join(r"test_objs", f"{fault.name}_depth_contours.obj")
    geojson_file = os.path.join(r"test_contours", f"{fault.name}_depth_contours.geojson")
    fault.contours.to_file(geojson_file, driver="GeoJSON")
    mesh.write(vtk_file)
    mesh.write(obj_file)



