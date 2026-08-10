
import pyvista as pv
from shapely.geometry import LineString
import numpy as np
from fault_mesh.faults.generic import calculate_strike_direction
import geopandas as gpd
import os

def contour_mesh(mesh_name: str, depth_interval: float = 500., 
                 z_range: tuple = (-2.e4, 0.), espg: int = 2193) -> gpd.GeoSeries:
    """
    Contour a mesh and save the contours as GeoJSON.

    Parameters:
    - mesh_name: Name of the mesh file.
    - depth_interval: Interval for contouring.
    - z_range: Range of Z values to consider for contours.
    - espg: EPSG code for the coordinate reference system.

    """
    assert os.path.exists(mesh_name), f"Mesh file {mesh_name} does not exist."
    mesh = pv.read(mesh_name)
    mesh['Z Height'] = mesh.points[:, 2]
    contours = mesh.contour(isosurfaces=int((z_range[1] - z_range[0]) / depth_interval) + 1, 
                            scalars='Z Height', rng=z_range)
    contour_levels = np.unique(contours.points[:, -1])

    contour_lines = []
    for level in contour_levels:
        contour = contours.points[np.isclose(contours.points[:, -1], level, atol=1e-3)]
        strike = calculate_strike_direction(contour[:, 0], contour[:, 1])
        strike_vector = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.])
        point_mean = contour.mean(axis=0)
        point_vectors = contour - point_mean
        distance_along_strike = np.dot(point_vectors, strike_vector)
        sorted_indices = np.argsort(distance_along_strike)
        sorted_contour = contour[sorted_indices]
        contour_line = LineString(sorted_contour)
        contour_lines.append(contour_line)

    contour_lines = gpd.GeoDataFrame({"depth": contour_levels}, 
                                     geometry=contour_lines, crs=espg)
    return contour_lines