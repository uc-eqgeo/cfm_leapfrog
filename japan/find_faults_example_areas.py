import geopandas as gpd
import pandas as pd
from fault_mesh.utilities.merging import merge_multiple_nearly_adjacent_segments

# Read in the fault traces
all_faults = gpd.read_file("japan_faults_utm_edited.geojson")
# Read in polygons that define areas of interest
example1_poly = gpd.read_file("example1.geojson")
example2_poly = gpd.read_file("example2.geojson")

# Trim the fault traces to the areas of interest
example1 = all_faults[all_faults.intersects(example1_poly.unary_union)]
example2 = all_faults[all_faults.intersects(example2_poly.unary_union)]

# Write the trimmed fault traces to file
example1.to_file("example1_faults.geojson", driver="GeoJSON")
example2.to_file("example2_faults.geojson", driver="GeoJSON")

# Read fault trace metadata from file
fault_metadata = pd.read_csv("segment_list_e.csv", encoding="unicode_escape", index_col=False)
# Trim metadata dataframe to the areas of interest
example1_metadata = fault_metadata[fault_metadata["number of behavioral segment"].isin(example1["Name"])]
example2_metadata = fault_metadata[fault_metadata["number of behavioral segment"].isin(example2["Name"])]

# Remove non-ascii characters from the dip metadata
dip_noascii = example1_metadata["dip"].str.encode("ascii", "ignore").str.decode("ascii")
dip_noascii = dip_noascii.str.split("  ", expand=True)
example1_metadata["dip"] = pd.to_numeric(dip_noascii.iloc[:, 0])
example1_metadata["dip_direction"] = dip_noascii.iloc[:, 1]

# Extract only columns of interest from the example 1 metadata dataframe
trimmed_example1 = example1_metadata[["number of behavioral segment", "name of behavioral segment", "dip", "dip_direction", "slip rate[m/ky]"]]
# Rename columns to make them easier to deal with
trimmed_example1.rename(columns={"number of behavioral segment": "fault_number", "name of behavioral segment": "fault_name", "slip rate[m/ky]": "slip_rate"}, inplace=True)

# Get fault trace geometries that match fault numbers in the metadata
example1_geom = all_faults[all_faults.Name.isin(trimmed_example1["fault_number"])]
# Get names and columns for faults
example1_geom = example1_geom[["Name", "geometry"]]
# Extract only long faults
example1_long = example1_geom[example1_geom.geometry.length > 7000.]

# Merge faults that have the same number but different trace geometries
included_faults = []
included_geoms = []
for fault in trimmed_example1["fault_number"]:
    matching = example1_long[example1_long.Name == fault]
    if len(matching) > 0:
        if len(matching.Name) > 1:
            merged_geom = merge_multiple_nearly_adjacent_segments(list(matching.geometry.explode()), tolerance=10000.)

        else:
            merged_geom = matching.geometry.explode().iloc[0]
        included_faults.append(fault)
        included_geoms.append(merged_geom)

# Extract fault metadata for faults that have been merged
output_df = trimmed_example1[trimmed_example1.fault_number.isin(included_faults)]
# Rename dataframe columns so they are equivalent to NZ CFM headings (can change later)
output_df.rename(columns={"fault_number": "Fault_ID", "fault_name": "Name", "dip": "Dip_pref", "dip_direction": "Dip_dir", "slip_rate": "SR_pref"}, inplace=True)
# Add geometries
output_gdf = gpd.GeoDataFrame(output_df, geometry=included_geoms, crs="EPSG:32654")
# Write to file (GeoJSON but could be shapefile)
output_gdf.to_file("example1_faults_trimmed.geojson", driver="GeoJSON")

# example1_metadata.to_csv("example1_metadata.csv", index=False)