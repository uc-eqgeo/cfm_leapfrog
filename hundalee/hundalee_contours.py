from fault_mesh.faults.leapfrog import LeapfrogFault
import geopandas as gpd
import numpy as np

for dip in [30., 40., 50., 60.]:
    fault = LeapfrogFault()
    trace = gpd.read_file("hundalee_trace.shp").geometry.values[0]
    fault.nztm_trace = trace

    fault.dip_dir_str = "NW"
    fault.dip_best = dip

    contour_intervals = np.arange(2000, 22000., 2000.)

    fault.generate_depth_contours(contour_intervals)
    fault.contours.to_file(f"contours/hundalee_contours{int(dip)}.shp")
    fault.calculate_footprint()
    gpd.GeoSeries(fault.footprint.boundary, crs=2193).to_file(f"footprints/hundalee_footprint{int(dip)}.shp")