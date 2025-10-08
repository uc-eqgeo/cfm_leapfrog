import pyvista as pv
from pathlib import Path
from shapely.geometry import LineString
import numpy as np
from fault_mesh.faults.generic import calculate_strike_direction
import geopandas as gpd
# import gemgis as gg

# fname = "./Alpine_full_segmented_transition_1_12_05km.stl"
# fname = "./no_dip_change_1km.stl"
# fname = "./truncated_south_1km.stl"
fname = "./Truncated_south.obj"
# end_points = gpd.read_file(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\howarth_alpine\endpoints\Wairau 3 combined_varying_dip_endpoints.shp")

mesh = pv.read(fname)
mesh['Z Height'] = mesh.points[:, 2]
contours = mesh.contour(isosurfaces=41, scalars='Z Height', rng=[-2.e4, 0.])
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

contour_lines = gpd.GeoSeries(contour_lines, crs="EPSG:2193")
contour_gdf = gpd.GeoDataFrame({"depth": contour_levels}, geometry=contour_lines, crs="EPSG:2193")
contour_gdf.to_file("alpine_contours_truncated_south.geojson", driver='GeoJSON')







pl = pv.Plotter()
pl.add_mesh(mesh, opacity=0.85)
pl.add_mesh(contours, color='white', line_width=5)
pl.show()