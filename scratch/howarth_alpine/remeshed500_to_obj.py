import meshio
import os
from pathlib import Path

# Create output directory if it doesn't exist
output_dir = Path('remeshed500_obj')
output_dir.mkdir(exist_ok=True)

# Get all STL files from input directory
input_dir = Path('remeshed500')
stl_files = list(input_dir.glob('*.stl'))

# Process each STL file
for stl_file in stl_files:
    # Read the STL file
    mesh = meshio.read(stl_file)
    
    # Create output filename with .obj extension
    output_file = output_dir / (stl_file.stem + '.obj')
    
    # Write as OBJ
    meshio.write(
        output_file,
        mesh,
        file_format="obj"
    )
    print(f"Converted {stl_file.name} to {output_file.name}")