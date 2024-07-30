import geopandas as gpd
import pandas as pd
from fault_mesh.utilities.merging import merge_multiple_nearly_adjacent_segments

# Read in the fault traces
# all_faults = gpd.read_file("japan_faults_utm_edited.geojson")
all_faults = gpd.read_file("faults_maehara_utm.geojson")

# Read fault trace metadata from file
fault_metadata = pd.read_csv("segment_list_e_wJSHISrev.csv", encoding="unicode_escape", index_col=False)
# # Trim metadata dataframe to the areas of interest
selected_metadata = fault_metadata[fault_metadata["number of behavioral segment"].isin(all_faults["Name"])].copy()
problem_metadata = selected_metadata[~selected_metadata["Dip_JSHIS"].notna()]
problem_metadata.to_csv("problem_metadata.csv", index=False)

selected_metadata = selected_metadata[selected_metadata["Dip_JSHIS"].notna()]

# Remove non-ascii characters from the dip metadata
# dip_noascii = selected_metadata["dip"].str.encode("ascii", "ignore").str.decode("ascii")
# dip_noascii = dip_noascii.str.split("  ", expand=True)
# selected_metadata["dip"] = pd.to_numeric(dip_noascii.iloc[:, 0])
# selected_metadata["dip_direction"] = dip_noascii.iloc[:, 1]

# Extract only columns of interest from the example 1 metadata dataframe
# trimmed_selected = selected_metadata[["number of behavioral segment", "name of behavioral segment", "dip", "dip_direction", "slip rate[m/ky]"]]
trimmed_selected = selected_metadata[["number of behavioral segment", "name of behavioral segment", "Dip_JSHIS", "DIP_JSHIS_dir", "slip rate[m/ky]"]]
# Rename columns to make them easier to deal with
# trimmed_selected = trimmed_selected.rename(columns={"number of behavioral segment": "fault_number", "name of behavioral segment": "fault_name", "slip rate[m/ky]": "slip_rate"})
trimmed_selected = trimmed_selected.rename(columns={"number of behavioral segment": "fault_number", "name of behavioral segment": "fault_name", "slip rate[m/ky]": "slip_rate",
                                                    "Dip_JSHIS": "dip", "DIP_JSHIS_dir": "dip_direction"})

# Get fault trace geometries that match fault numbers in the metadata
selected_geom = all_faults[all_faults.Name.isin(trimmed_selected["fault_number"])]
# Get names and columns for faults
selected_geom = selected_geom[["Name", "geometry"]]
# Extract only long faults
selected_long = selected_geom[selected_geom.geometry.length > 7000.]

# Merge faults that have the same number but different trace geometries
included_faults = []
included_geoms = []
for fault in trimmed_selected["fault_number"]:
    matching = selected_long[selected_long.Name == fault]
    if len(matching) > 0:
        if len(matching.Name) > 1:
            merged_geom = merge_multiple_nearly_adjacent_segments(list(matching.geometry.explode(index_parts=True)), tolerance=20000.)

        else:
            merged_geom = matching.geometry.explode(index_parts=True).iloc[0]
        included_faults.append(fault)
        included_geoms.append(merged_geom)

# Extract fault metadata for faults that have been merged
output_df = trimmed_selected[trimmed_selected.fault_number.isin(included_faults)]
# Rename dataframe columns so they are equivalent to NZ CFM headings (can change later)
output_df = output_df.rename(columns={"fault_number": "Fault_ID", "fault_name": "Name", "dip": "Dip_pref", "dip_direction": "Dip_dir", "slip_rate": "SR_pref"})
# Add geometries
output_gdf = gpd.GeoDataFrame(output_df, geometry=included_geoms, crs="EPSG:32654")
# Write to file (GeoJSON but could be shapefile)
output_gdf.to_file("selected_faults_trimmed.geojson", driver="GeoJSON")

# selected_metadata.to_csv("selected_metadata.csv", index=False)