from rsqsim_api.fault.multifault import RsqSimSegment
from rsqsim_api.fault.splines import fault_edge_spline, fault_depth_spline
from rsqsim_api.io.array_operations import read_gmt_grid
import numpy as np
from scipy.spatial import KDTree
import meshio
from glob import glob
import os.path

# Parameters for smoothing
fault_slip_spline = fault_edge_spline(1., 7500., 300000.,
                                        min_slip_fraction=0., line_stop_fraction=0.5, gradient_change=1.5,
                                        spline_k=3, spline_s=0.0, resolution=10.)
depth_slip_spline = fault_depth_spline(0.6, 0.8, after_change_gradient=1.15)

depth_x, depth_y, depth_z = read_gmt_grid("with_hannu_mods_resampled.grd")
depth_X, depth_Y = np.meshgrid(depth_x, depth_y)
depth_tree = KDTree(np.column_stack([depth_X.ravel(), depth_Y.ravel(), depth_z.ravel()]))


in_files = list(glob("./remeshed_faults/*.vtk"))
out_dir = "./remeshed_tapered_slip"
for file in in_files:
    basename = os.path.basename(file).split(".")[0]
    print(basename)
    fault = RsqSimSegment.from_vtk(file)
    # fault.to_vtk(os.path.join(out_dir, basename + "_collapsed.vtk"), write_slip=True)
    fault_bottom_edge = fault.bottom_edge_point_cloud(num_interp=25, depth_tolerance=100.)
    fault_bottom_edge_tree = KDTree(fault_bottom_edge)
    fault_centres = fault.get_patch_centres()
    distances, _ = fault_bottom_edge_tree.query(fault_centres)
    _, depth_indices = depth_tree.query(fault_centres)
    depth_values = depth_z.ravel()[depth_indices]
    depth_proportions = fault_centres[:, 2] / depth_values
    depth_proportions[depth_proportions > 1.] = 1.
    depth_proportions[depth_proportions < 0.] = 0.
    depth_multiplier = depth_slip_spline(depth_proportions)
    slip_multiplier = fault_slip_spline(distances)
    min_slip_multiplier = np.minimum(depth_multiplier, slip_multiplier)

    mesh = meshio.read(file)
    mesh_slip_rate = mesh.cell_data["slip"][0]
    new_slip_rate = mesh_slip_rate * min_slip_multiplier
    mesh.cell_data["slip"][0] = new_slip_rate
    meshio.write(os.path.join(out_dir, basename + "_tapered_slip.vtk"), mesh)


