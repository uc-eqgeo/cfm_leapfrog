from eq_fault_geom.geomio.cfm_faults import CfmMultiFault
import geopandas as gpd

shp = "../gis/cfm_gt_1_5.gpkg"
gdf = gpd.read_file(shp)
data_d90 = CfmMultiFault.from_shp(shp, depth_type="D90", sort_sr=True)

f0 = data_d90.faults[0]
