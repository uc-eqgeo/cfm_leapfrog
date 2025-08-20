from fault_mesh.utilities.meshing import get_strike_dip_from_normal, fit_plane_to_points
from fault_mesh.utilities.splines import spline_fit_contours
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString
from scipy.interpolate import RBFInterpolator 
from matplotlib import pyplot as plt
import triangle as tr
import meshio

# contours = gpd.read_file("C:\\Users\\arh128\\PycharmProjects\\cfm_leapfrog\\scratch\\cfm_no_leapfrog\\NZ_CFM_v1_0_rs1km_modified_GDM_Final\\depth_contours\\Alpine George - Wairau combined_contours.shp")
contours = gpd.read_file("C:\\Users\\arh128\\PycharmProjects\\cfm_leapfrog\\scratch\\cfm_no_leapfrog\\NZ_CFM_v1_0_rs1km_modified_GDM_Final\\depth_contours\\Wellington Hutt Valley combined_contours.shp")
spline_contours = spline_fit_contours(contours, point_spacing=100., output_spacing=1000.)
segmentized = spline_contours.segmentize(1000.)
spline_contour_list = [np.array(contour.coords) for contour in segmentized.geometry.values]
spline_contour_points = np.vstack(spline_contour_list)
spline_contour_dict = {tuple(point): i for (i,point) in enumerate(spline_contour_points)}

spline_contour_boundary = np.vstack([spline_contour_list[0], 
                            np.vstack([spline_contour_list[i][-1] for i in range(1, len(spline_contour_list) -1)]),
                            spline_contour_list[-1][::-1],
                            np.vstack([spline_contour_list[i][0] for i in range(len(spline_contour_list) -1, 0, -1)])])

contours_list = []
for contour in contours.geometry.values:
    if isinstance(contour, LineString):
        contours_list.append(np.array(contour.coords))
    else:
        assert isinstance(contour, MultiLineString)
        contours_list.append(np.vstack([np.array(line.coords) for line in contour.geoms]))


all_contour_points = np.vstack(contours_list)
all_contour_dict = {tuple(point): i for (i,point) in enumerate(all_contour_points)}
plane_normal, plane_origin = fit_plane_to_points(all_contour_points, eps=1.0e-5)

# Calculate strike and dip from normal
strike, dip = get_strike_dip_from_normal(plane_normal)
print(f"Plane strike: {strike:.1f}°, dip: {dip:.1f}°")

# Calculate in-plane x vector (along strike)
strike_rad = np.radians(strike)
in_plane_x_vector = np.array([np.sin(strike_rad), np.cos(strike_rad), 0])
in_plane_x_vector = in_plane_x_vector - np.dot(in_plane_x_vector, plane_normal) * plane_normal
in_plane_x_vector = in_plane_x_vector / np.linalg.norm(in_plane_x_vector)

# Calculate in-plane y vector (down dip)
in_plane_y_vector = np.cross(plane_normal, in_plane_x_vector)
in_plane_y_vector = in_plane_y_vector / np.linalg.norm(in_plane_y_vector)


# resolve each contour in the list into the plane
resolved_contours = []
for contour in contours_list:
    # Translate the contour to the origin
    contour -= plane_origin
    # Rotate the contour to align with the in-plane x vector
    rotated_contour_x = np.dot(contour, in_plane_x_vector)
    rotated_contour_y = np.dot(contour, in_plane_y_vector)
    rotated_contour_z = np.dot(contour, plane_normal)

    # Create a new contour with the rotated coordinates
    resolved_contours.append(np.vstack([rotated_contour_x, rotated_contour_y, rotated_contour_z]).T)

all_resolved_points = np.vstack(resolved_contours)
np.savetxt("all_resolved_points.txt", all_resolved_points)

# resolve the contour boundary into the plane
resolved_boundary_x = np.dot(spline_contour_boundary, in_plane_x_vector)
resolved_boundary_y = np.dot(spline_contour_boundary, in_plane_y_vector)
resolved_boundary = np.vstack([resolved_boundary_x, resolved_boundary_y]).T

resolved_spline_contours = []
for contour in spline_contour_list:
    # Translate the contour to the origin
    contour -= plane_origin

    # Rotate the contour to align with the in-plane x vector
    rotated_contour_x = np.dot(contour, in_plane_x_vector)
    rotated_contour_y = np.dot(contour, in_plane_y_vector)

    # Create a new contour with the rotated coordinates
    resolved_spline_contours.append(np.vstack([rotated_contour_x, rotated_contour_y]).T)

all_resolved_spline_points = np.vstack(resolved_spline_contours)

interpolator = RBFInterpolator(all_resolved_points[:, :2], all_resolved_points[:, 2], kernel='thin_plate_spline')

interpolated = interpolator(all_resolved_spline_points)

# Turn back into contour coordinates
interpolated_xyz = plane_origin + all_resolved_spline_points[:, 0][:, np.newaxis] * np.repeat([in_plane_x_vector], interpolated.shape[0], axis=0) + all_resolved_spline_points[:, 1][:, np.newaxis] * np.repeat([in_plane_y_vector], interpolated.shape[0], axis=0) + interpolated[:, np.newaxis] * np.repeat([plane_normal], interpolated.shape[0], axis=0)

plt.close("all")
fig, ax = plt.subplots()
ax.plot(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", label="Resolved boundary")
ax.scatter(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", s=1, label="Resolved boundary points")
plt.show(block=True)

boundary_segments = np.array([[spline_contour_dict[tuple(spline_contour_boundary[i])], spline_contour_dict[tuple(spline_contour_boundary[i + 1])]] for i in range(len(spline_contour_boundary) - 1)])
boundary_segments = np.vstack([boundary_segments, [spline_contour_dict[tuple(spline_contour_boundary[-1])], spline_contour_dict[tuple(spline_contour_boundary[0])]]])
A = dict(vertices=all_resolved_spline_points, segments=boundary_segments)
B = tr.triangulate(A, 'p')

plt.close("all")
tr.compare(plt, A, B)
plt.show(block=True)

# write the mesh to a file
mesh = meshio.Mesh(points=interpolated_xyz, cells={"triangle": B['triangles']})
mesh.write("alpine_rbf.vtk", file_format="vtk")

