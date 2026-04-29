import pyvista as pv
import meshio

from glob import glob
meshes = [pv.from_meshio(meshio.read(f)) for f in glob(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\revised_cutting\vtks\*.vtk")]
combined = pv.merge(meshes)
combined.save(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\merged_mesh_uncut.vtk")

for f in glob(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\revised_cutting\vtks\*.vtk"):
    stl_name = f.replace(".vtk", ".stl").replace("vtks", "stls")
    meshio.write(stl_name, meshio.read(f), file_format="stl")

meshes = [pv.from_meshio(meshio.read(f)) for f in glob(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\revised_cutting\cut_meshes\*.vtk")]
combined = pv.merge(meshes)
combined.save(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\merged_mesh_cut.vtk")
