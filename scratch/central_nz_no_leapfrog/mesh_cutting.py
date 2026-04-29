from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.mesh import FaultMesh
from fault_mesh.io.array_operations import read_raster
import pyvista as pv
import os

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

# Read depth raster and convert to PyVista surface
depth_pyvista = read_raster("with_hannu_mods.tif", use_z=True, out_crs=f"EPSG:{epsg}")

fault_data.read_additional_cuts("additional_cuts.csv")

# Read in existing meshes for each fault if they exist, and store in a dictionary for cutting
for fault in fault_data.curated_faults:
    obj_file = os.path.join(r"test_objs", f"{fault.name}_depth_contours.obj")
    if os.path.exists(obj_file):
        print(f"Reading existing mesh for {fault.name}")
        mesh = FaultMesh.from_file(obj_file)
        fault.mesh = mesh

# Dictionary for use with cutting hierarchy
cutting_dict = {fault.name: fault for fault in fault_data.curated_faults if fault.mesh is not None}

# Create output directory for final cut meshes
output_dir = r"test_final_meshes"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Threshold currently does nothing, ignore for now.
threshold = 0.7
# Set a minimum distance threshold to avoid cutting faults that are very close together (e.g. due to KDTree artefacts or small faults right next to each other)
min_distance = 2.e3
for fault_name in fault_data.cutting_hierarchy:
    print (f"Processing {fault_name}...")
    fault = fault_data.name_dict[fault_name]
    # Get list of faults that cut this fault, in order of hierarchy
    faults_to_cut =  fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(fault_name)]
    # Loop through cutting faults and decide whether to cut based on intersection with cutting fault mesh and distance to cutting fault trace. 
    for cut_name in faults_to_cut:
        if cut_name in cutting_dict:
            higher_faults_to_cut = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(cut_name)]
            higher_meshes = [cutting_dict[name].mesh for name in higher_faults_to_cut]
            cutting_mesh = cutting_dict[cut_name].mesh
            decide_to_cut = fault.mesh.decide_whether_to_cut(cutting_mesh, threshold=threshold, min_distance=min_distance, 
                                                             higher_meshes=higher_meshes, bottom_depth=-31500., fancy_cutting=True)
            if decide_to_cut:
                print(f"Cutting {fault.name} by {cut_name}")
                cutting_fragment = cutting_mesh.generate_cutting_mesh(fault.mesh, max_distance=5.e3)
                new_mesh = fault.mesh.cut_mesh(cutting_mesh, fault_trace=fault.original_nztm_trace_array, cutting_fragment=cutting_fragment, fancy_cutting=True)
                fault_data.name_dict[fault.name].mesh = new_mesh
                cutting_dict[fault.name].mesh = new_mesh

    # Cut meshes by depth surface and save final meshes as OBJ files     
    try:
        depth_trimmed = fault_data.name_dict[fault.name].mesh.cut_mesh_pv(depth_pyvista, fault_trace=fault.original_nztm_trace_array)
        fault_data.name_dict[fault.name].mesh = depth_trimmed
    except Exception as e:
        print(f"Error trimming depth for {fault.name}: {e}")
        continue
    output_file = os.path.join(output_dir, f"{fault.name}_cut.obj")
    fault_data.name_dict[fault.name].mesh.mesh.write(output_file)
    print(f"Done: {fault.name}")