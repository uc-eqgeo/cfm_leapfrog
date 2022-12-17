import geopandas as gpd
from matplotlib import pyplot as plt
cfm = gpd.read_file("cfm_central_nz.geojson")
cfm_gt1_5 = cfm[cfm.SR_pref >= 1.5]
cfm_sr = cfm_gt1_5.sort_values(by="SR_pref", ascending=False, ignore_index=True)
print(cfm_sr.Name)

cfm_sr.to_file("cfm_gt_1_5.geojson", driver="GeoJSON")

