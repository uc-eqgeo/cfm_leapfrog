import pyvista as pv
import meshio

from glob import glob
meshes = [pv.from_meshio(meshio.read(f)) for f in glob(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\*.vtk")]
combined = pv.merge(meshes)
combined.save(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\merged_mesh.vtk")