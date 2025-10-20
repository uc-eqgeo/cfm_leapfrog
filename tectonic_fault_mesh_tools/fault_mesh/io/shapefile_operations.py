import geopandas as gpd
from typing import Union
from shapely.geometry import Polygon, LineString, MultiLineString, MultiPolygon, Point, MultiPoint
import os
import pandas as pd


def geopandas_polygon_to_gmt(gdf: Union[gpd.GeoDataFrame, gpd.GeoSeries], out_file: str, matlab: bool = False):
    """
    :param gdf: Geodataframe with polygon or multipolygon geometry
    :param out_file: normally has a ".gmt" ending
    :return:
    """
    # Check that writing polygon
    assert all([type(x) in (Polygon, MultiPolygon) for x in gdf.geometry])
    out_id = open(out_file, "w")
    geom_ls = []
    for geom in gdf.geometry:
        if isinstance(geom, Polygon):
            geom_ls.append(geom)
        elif isinstance(geom, MultiPolygon):
            geom_ls += list(geom)
    for geom in geom_ls:
        if isinstance(geom.boundary, LineString):
            x, y = geom.boundary.coords.xy
            if matlab:
                out_id.write("NaN\n")
            else:
                out_id.write(">\n")
            for xi, yi in zip(x, y):
                out_id.write("{:.4f} {:.4f}\n".format(xi, yi))
        elif isinstance(geom.boundary, MultiLineString):
            for line in geom.boundary:
                x, y = line.coords.xy
                if matlab:
                    out_id.write("NaN\n")
                else:
                    out_id.write(">\n")
                for xi, yi in zip(x, y):
                    out_id.write("{:.4f} {:.4f}\n".format(xi, yi))
    out_id.close()
    return


def geopandas_linestring_to_gmt(gdf: Union[gpd.GeoDataFrame, gpd.GeoSeries], out_file: str, matlab: bool = False):
    """

    :param gdf: Geodataframe or geoseries with polygon or multipolygon geometry
    :param out_file: normally has a ".gmt" ending
    :return:
    """
    # Check that writing lines
    assert all([type(x) in (LineString, MultiLineString) for x in gdf.geometry])
    out_id = open(out_file, "w")
    geom_ls = []
    for geom in gdf.geometry:
        if isinstance(geom, LineString):
            geom_ls.append(geom)
        elif isinstance(geom, MultiLineString):
            geom_ls += list(geom)
    for geom in geom_ls:
        if isinstance(geom, LineString):
            x, y = geom.coords.xy
            if matlab:
                out_id.write("NaN\n")
            else:
                out_id.write(">\n")
            for xi, yi in zip(x, y):
                out_id.write("{:.4f} {:.4f}\n".format(xi, yi))
        elif isinstance(geom, MultiLineString):
            for line in geom:
                x, y = line.coords.xy
                if matlab:
                    out_id.write("NaN\n")
                else:
                    out_id.write(">\n")
                for xi, yi in zip(x, y):
                    out_id.write("{:.4f} {:.4f}\n".format(xi, yi))
    out_id.close()
    return


def geopandas_points_to_gmt(gdf: Union[gpd.GeoDataFrame, gpd.GeoSeries], out_file: str):
    """

    :param gdf: Geodataframe or geoseries with point or multipoint geometry
    :param out_file: normally has a ".gmt" ending
    :return:
    """
    # Check that writing points
    assert all([type(x) in (Point, MultiPoint) for x in gdf.geometry])
    out_id = open(out_file, "w")
    geom_ls = []
    for geom in gdf.geometry:
        if isinstance(geom, Point):
            geom_ls.append(geom)
        elif isinstance(geom, MultiPoint):
            geom_ls += list(geom)
    for geom in geom_ls:
        if isinstance(geom, Point):
            out_id.write("{:.4f} {:.4f}\n".format(geom.x, geom.y))
        elif isinstance(geom, MultiPoint):
            for x, y in geom.xy:
                out_id.write("{:.4f} {:.4f}\n".format(x, y))
    out_id.close()
    return


def shp_to_gmt(shapefile: str, out_file: str, in_epsg: int = None, out_epsg: int = None):
    """

    :param shapefile:
    :param out_file:
    :param in_epsg:
    :param out_epsg:
    :return:
    """
    assert os.path.exists(shapefile)
    if in_epsg is not None:
        gdf = gpd.GeoDataFrame.from_file(shapefile, crs={"init": "{:d}".format(in_epsg)})
    else:
        gdf = gpd.GeoDataFrame.from_file(shapefile)

    if out_epsg is not None:
        if in_epsg != out_epsg:
            out_gdf = gdf.to_crs(epsg=out_epsg)
        else:
            out_gdf = gdf

    else:
        out_gdf = gdf

    if any([isinstance(gdf.geometry[0], x) for x in (Point, MultiPoint)]):
        geopandas_points_to_gmt(out_gdf, out_file)
    elif any([isinstance(gdf.geometry[0], x) for x in (LineString, MultiLineString)]):
        geopandas_linestring_to_gmt(out_gdf, out_file)
    elif any([isinstance(gdf.geometry[0], x) for x in (Polygon, MultiPolygon)]):
        geopandas_polygon_to_gmt(out_gdf, out_file)
    else:
        raise TypeError("shapefile must contain (only) one of points, lines or polygons")

    return


def reproject_csv(in_csv: str, out_csv: str = None, sep: str = None, in_epsg: int = 2193, out_epsg: int = 4326,
                  header: int = None):
    """

    :param in_csv:
    :param out_csv:
    :param sep:
    :param in_epsg:
    :param out_epsg:
    :param header:
    :return:
    """
    assert in_epsg != out_epsg
    assert os.path.exists(in_csv)
    if sep is None:
        df = pd.read_csv(in_csv, delim_whitespace=True, header=header)
    else:
        df = pd.read_csv(in_csv, sep=sep, header=header)

    assert len(df.iloc[0]) > 2, "Expecting at least 2 columns (lon, lat)"

    geom = gpd.GeoSeries([Point(xi, yi) for xi, yi in zip(df.iloc[:, 0], df.iloc[:, 1])],
                         crs={"init": "epsg:{:d}".format(in_epsg)})
    new_geom = geom.to_crs(epsg=out_epsg)

    df.iloc[:, 0] = new_geom.x
    df.iloc[:, 1] = new_geom.y

    if out_csv is None:
        if "." in in_csv:
            end = "." + in_csv.split(".")[-1]
            start = in_csv.split(end)[0]
        else:
            start = in_csv
        if out_epsg == 2193:
            out_name = start + "_nztm.csv"
        elif out_epsg == 4326:
            out_name = start + "_wgs.csv"
        else:
            out_name = start + "_reprojected.csv"
    else:
        out_name = out_csv

    out_sep = sep if sep is not None else " "
    df.to_csv(out_name, sep=out_sep, header=header, index=None)



