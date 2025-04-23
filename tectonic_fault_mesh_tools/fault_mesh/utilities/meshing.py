"""
Utilities for meshing from fault contours.
This module contains functions for triangulating contours and generating meshes from them.
Functions:
    - get_strike_dip_from_normal: Calculate strike and dip from a normal vector
    - fit_plane_to_points: Fit a plane to a set of 3D points using SVD
    - triangulate_contours: Triangulate contours and write to a file
"""
import numpy as np
import matplotlib.pyplot as plt
import triangle as tr
import geopandas as gpd
import meshio

def get_strike_dip_from_normal(normal_vector):
    """
    :param normal_vector: The normal vector to the plane [nx, ny, nz]
    :type normal_vector: numpy.ndarray
    :return: A tuple containing (strike, dip) in degrees. Strike is measured 0-360 degrees clockwise from North. 
             Dip is measured 0-90 degrees from horizontal.
    :rtype: tuple(float, float)
    :note: Uses the right-hand rule convention for strike/dip measurements
    :example:
        >>> normal = np.array([0, 0, 1])
        >>> strike, dip = get_strike_dip_from_normal(normal) 
        >>> print(strike, dip)
        0.0, 0.0
    """
    nx, ny, nz = normal_vector
    
    # Calculate dip (angle from horizontal)
    dip = np.degrees(np.arccos(abs(nz)))
    
    # Calculate strike (direction of horizontal line in plane)
    # Strike is perpendicular to the horizontal component of the normal vector
    if nz == 0:  # Vertical plane
        strike = np.degrees(np.arctan2(-nx, ny))
    else:
        strike = np.degrees(np.arctan2(nx, ny)) - 90
    
    # Adjust strike to 0-360 range
    while strike < 0:
        strike += 360
    while strike >= 360:
        strike -= 360
    
    return strike, dip

def fit_plane_to_points(points: np.ndarray, eps: float=1.0e-5):
    """
    Find best-fit plane through a set of points using SVD.
    
    This function fits a plane to a set of 3D points by computing the singular value decomposition (SVD)
    of the moment matrix. The plane passes through the centroid of all points, which improves stability.
    
    :param points: Array of 3D points with shape (n, 3) where n is the number of points.
    :type points: numpy.ndarray
    :param eps: Threshold to zero out very small components of the normal vector, default is 1.0e-5.
    :type eps: float
    
    :return: A tuple containing (plane_normal, plane_origin)
             plane_normal is the unit normal vector to the fitted plane (A, B, C)
             plane_origin is the centroid of the input points, used as the plane origin
    :rtype: tuple(numpy.ndarray, numpy.ndarray)
    
    :note: The function uses SVD on a 3x3 covariance matrix rather than on the full point array,
           which is more computationally efficient. The returned normal vector is normalized and
           oriented so that the z-component is positive (when not zero).
    """
    # Compute plane origin and subract it from the points array.
    plane_origin = np.mean(points, axis=0)
    x = points - plane_origin

    # Dot product to yield a 3x3 array.
    moment = np.dot(x.T, x)

    # Extract single values from SVD computation to get normal.
    plane_normal = np.linalg.svd(moment)[0][:,-1]
    small = np.where(np.abs(plane_normal) < eps)
    plane_normal[small] = 0.0
    plane_normal /= np.linalg.norm(plane_normal)
    if (plane_normal[-1] < 0.0):
        plane_normal *= -1.0

    return plane_normal, plane_origin

def triangulate_contours(contours: gpd.GeoDataFrame, mesh_name: str, mesh_format="vtk", check_mesh: bool=False, check_strike_dip: bool=True):
    """
    Triangulate the contours and write to a file
    :param contours:
        The contours to triangulate
    :param mesh_name:                       
        The name of the mesh file to write to
    :param mesh_format:
        The format of the mesh file to write to
    :param check_mesh:
        Whether to check the mesh and plot it
    :param check_strike_dip:
        Whether to check the strike and dip of the plane
    :return:
        None
    """
    contours_list = [np.array(contour.coords) for contour in contours.geometry.values]
    contour_boundary = np.vstack([contours_list[0], 
                                np.vstack([contours_list[i][-1] for i in range(1, len(contours_list) -1)]),
                                contours_list[-1][::-1],
                                np.vstack([contours_list[i][0] for i in range(len(contours_list) -1, 0, -1)])])
    all_contour_points = np.vstack(contours_list)
    all_contour_dict = {tuple(point): i for (i,point) in enumerate(all_contour_points)}
    plane_normal, plane_origin = fit_plane_to_points(all_contour_points, eps=1.0e-5)

    # Calculate strike and dip from normal
    strike, dip = get_strike_dip_from_normal(plane_normal)
    if check_strike_dip:
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
        # Create a new contour with the rotated coordinates
        resolved_contours.append(np.vstack([rotated_contour_x, rotated_contour_y]).T)

    all_resolved_points = np.vstack(resolved_contours)

    # resolve the contour boundary into the plane
    resolved_boundary_x = np.dot(contour_boundary, in_plane_x_vector)
    resolved_boundary_y = np.dot(contour_boundary, in_plane_y_vector)
    resolved_boundary = np.vstack([resolved_boundary_x, resolved_boundary_y]).T

    if check_mesh:
        plt.close("all")
        fig, ax = plt.subplots()
        ax.plot(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", label="Resolved boundary")
        ax.scatter(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", s=1, label="Resolved boundary points")
        plt.show(block=True)
    # mesh the resolved contours
    boundary_segments = np.array([[all_contour_dict[tuple(contour_boundary[i])], all_contour_dict[tuple(contour_boundary[i + 1])]] for i in range(len(contour_boundary) - 1)])
    boundary_segments = np.vstack([boundary_segments, [all_contour_dict[tuple(contour_boundary[-1])], all_contour_dict[tuple(contour_boundary[0])]]])
    A = dict(vertices=all_resolved_points, segments=boundary_segments)
    B = tr.triangulate(A, 'p')

    if check_mesh:
        plt.close("all")
        tr.compare(plt, A, B)
        plt.show(block=True)

    # write the mesh to a file
    mesh = meshio.Mesh(points=all_contour_points, cells={"triangle": B['triangles']})
    mesh.write(mesh_name, file_format=mesh_format)