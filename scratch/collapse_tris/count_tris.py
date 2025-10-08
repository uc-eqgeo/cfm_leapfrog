import meshio
import glob

stls = list(glob.glob("all_collapsed2500/*.stl"))
tris = 0
for stl in stls:
    mesh = meshio.read(stl, file_format="stl")
    tris += mesh.cells_dict["triangle"].shape[0]