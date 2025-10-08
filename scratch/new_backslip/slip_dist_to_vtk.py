import meshio
import numpy as np
from fault_mesh.utilities.splines import fault_depth_spline
from scipy.spatial import KDTree

# mesh = meshio.read("variable_dip_1km.stl")
mesh = meshio.read("truncated_south_1km.stl")
points = mesh.points
triangles = mesh.cells_dict["triangle"]
tri_verts = points[triangles]
centroids = tri_verts.mean(axis=1)
depths = centroids[:, 2]

data_for_interp = np.load("matched_slip_contours_truncated_south.npy")

depth_slip_spline = fault_depth_spline(0.6, 0.8, after_change_gradient=1.15)

data_points = data_for_interp[:, :3]
slip_rate = data_for_interp[:, 3]
rake = data_for_interp[:, 4]

data_tree = KDTree(data_points)

centroid_indices = data_tree.query(centroids)[1]
slip_rate_values = slip_rate[centroid_indices]
rake_values = rake[centroid_indices]
rake_values[rake_values > 360.] = -110.  # Ensure no negative rake values
rake_values[np.abs(rake_values - -180.) < 1e-6] = 180.  # Correct -180 to 180

slip_rate_with_depth = slip_rate_values * depth_slip_spline(depths / -20000.)

mesh.cell_data["slip"] = [slip_rate_with_depth]
mesh.cell_data["rake"] = [rake_values]

mesh.write("slip_distribution_truncated_south.vtk")