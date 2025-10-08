import geopandas as gpd
from pathlib import Path

from fault_mesh.utilities.meshing import triangulate_contours
from fault_mesh.utilities.splines import spline_fit_contours





# hope_combined = Path(__file__).parent / "depth_contours" / "Hope combined_contours.shp"
# hope_combined = gpd.read_file(hope_combined)

# spline_contours = spline_fit_contours(hope_combined, point_spacing=100., output_spacing=1000.)

# file_name = Path(__file__).parent / "depth_contours" / "resolved_contours.vtk"
# triangulate_contours(spline_contours, file_name, mesh_format="vtk", check_mesh=True)

alpine_combined = Path(__file__).parent / "depth_contours" / "Alpine Kaniere to Springs Junction combined_contours.shp"
alpine_combined = gpd.read_file(alpine_combined)

spline_contours = spline_fit_contours(alpine_combined, point_spacing=100., output_spacing=1000.)
file_name = Path(__file__).parent / "depth_contours" / "resolved_contours_alpine.vtk"
triangulate_contours(spline_contours, file_name, mesh_format="vtk", check_mesh=True)





