import os
import meshio
import glob

# Find all remeshed STL files in the current directory
stl_files = glob.glob("all_collapsed2500/*remeshed_collapsed.stl")

# Convert each STL file to OBJ
for stl_file in stl_files:
    # Read the STL file
    mesh = meshio.read(stl_file)
    
    # Create output filename by replacing .stl with .obj
    obj_file = "all_remeshed_obj/" + os.path.splitext(os.path.basename(stl_file))[0] + ".obj"
    
    # Write to OBJ format
    meshio.write(obj_file, mesh, file_format="obj")
    print(f"Converted {stl_file} to {obj_file}")