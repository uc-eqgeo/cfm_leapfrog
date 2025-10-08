import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import meshio

from fault_mesh.faults.mesh import FaultMesh

fname = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\Whitemans Valley_depth_contours.vtk"
alpine_fname = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\Alpine George - Wairau combined_depth_contours.vtk"
whitemans = pv.read(fname).extract_surface()
wmesh = meshio.read(fname)

fmesh = FaultMesh.from_file(alpine_fname)



# aotea = pv.read(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\Aotea - Evans Bay_depth_contours.vtk").extract_surface()

# alpine = pv.read(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\cfm_no_leapfrog\test_vtks\Alpine George - Wairau combined_depth_contours.vtk").extract_surface()

# collision, _ = whitemans.collision(aotea, generate_scalars=True)
# scalars = np.zeros(collision.n_cells, dtype=bool)
# scalars[collision.field_data['ContactCells']] = True

# clipped = aotea.clip_surface(whitemans, invert=False)

# p = pv.Plotter()
# # p.add_mesh(whitemans, color='red', opacity=0.5, show_edges=True)
# p.add_mesh(aotea, color='blue', opacity=0.5, show_edges=True)
# p.add_mesh(collision, scalars=scalars)
# # p.add_mesh(clipped, color='green', opacity=0.5, show_edges=True)

# p.show()