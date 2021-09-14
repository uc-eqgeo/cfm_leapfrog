import geopandas as gpd
import numpy as np
from shapely.ops import linemerge, unary_union, polygonize

rect = gpd.read_file("clipping_area.gpkg")
rect_poly = list(rect.explode().geometry)[0]

lines = list(gpd.read_file("clipping_lines.gpkg").geometry)
depth_points = gpd.read_file("depth_points.gpkg")

merged = linemerge([rect_poly.boundary] + lines)
union = unary_union(merged)
polys = list(polygonize(union))

depth_list = []
for poly in polys:
    for i, row in depth_points.iterrows():
        if row.geometry.within(poly):
            depth_list.append(row.depth)

polys_with_depths = gpd.GeoDataFrame({"depth": depth_list}, geometry=polys, crs=2193)
polys_with_depths.to_file("depth_polys_with_depths.gpkg", driver="GPKG")

with open("depth_polygons.gmt", "w") as fid:
    for i, row in polys_with_depths.iterrows():
        if row.depth > 0:
            fid.write(f"> -Z{row.depth * -1000}\n")
            for coord in np.array(row.geometry.boundary):
                fid.write(f"{coord[0]:.4f} {coord[1]:.4f}\n")

