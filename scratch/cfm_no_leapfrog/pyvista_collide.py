import pyvista as pv
from pathlib import Path
from fault_mesh.io.array_operations import read_raster

# read in DEM
dem_name = r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\depth_and_interface\ellis_depth_surface\with_hannu_mods.tif"
dem = read_raster(dem_name, out_crs="EPSG:2193", use_z=True)
dem_surface = dem.extract_surface()

hope_name = Path(__file__).parent / "depth_contours" / "resolved_contours.vtk"
alpine_name = Path(__file__).parent / "depth_contours" / "resolved_contours_alpine.vtk"
hope_mesh = pv.read(hope_name).extract_surface()
alpine_mesh = pv.read(alpine_name).extract_surface()

# results = alpine_mesh.intersection(hope_mesh)

# clipped_mesh = hope_mesh.clip_surface(alpine_mesh, invert=False)
# clipped_mesh = clipped_mesh.clip_surface(dem_surface, invert=False)

distances = hope_mesh.compute_implicit_distance(alpine_mesh)
# inner = distances.threshold((-1.0, 1.0), scalars="implicit_distance", continuous=True)
# inner = distances.extract_values(0., scalars="implicit_distance", preference="cell", include_cells=True, component_mode="any")
pl = pv.Plotter()
pl.add_mesh(alpine_mesh, color="blue", opacity=0.5, style="wireframe")
edges = alpine_mesh.extract_feature_edges()
pl.add_mesh(edges, color="red", opacity=0.5, style="wireframe")
# pl.add_mesh(inner, opacity=0.5, clim=[-100, 100])
# pl.add_mesh(clipped_mesh, color="red", opacity=0.5, style="wireframe")
# pl.add_mesh(dem_surface, color="white", opacity=0.5, style="wireframe")
pl.show()