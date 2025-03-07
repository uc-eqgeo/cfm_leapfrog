from fault_mesh.faults.leapfrog import LeapfrogMultiFault
import numpy as np

# Load the model:
model = LeapfrogMultiFault.from_shp(r"F:\vscode_projects\synscarp\taieri\titri_erin_trace_500m.shp", epsg=2193)
titri = model.faults[0]

# Generate the contours:
contour_depths = np.arange(0, 28000, 2000)
contours = titri.generate_depth_contours(contour_depths, km=False)
titri.contours.to_file(r"F:\vscode_projects\synscarp\taieri\titri_contours_dip60.shp")

titri.footprint_geoseries.boundary.to_file(r"F:\vscode_projects\synscarp\taieri\titri_footprint_dip60.shp")