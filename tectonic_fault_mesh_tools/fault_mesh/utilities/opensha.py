from typing import Union
from shapely.geometry import LineString, Polygon
from xml.etree import ElementTree as ElemTree


def fault_trace_xml(geometry: LineString, section_name: str, z: Union[float, int] = 0):
    """
    XML element
    :param geometry: should be in lon lat
    :param section_name:
    :param z: Generally zero
    :return:
    """
    trace_element = ElemTree.Element("FaultTrace", attrib={"name": section_name})
    ll_float_str = "{:.4f}"
    # extract arrays of lon and lat
    x, y = geometry.xy
    # Loop through addis each coordinate as sub element
    for x_i, y_i in zip(x, y):
        if x_i <= 0.:
            x_i += 360.
        loc_element = ElemTree.Element("Location", attrib={"Latitude": ll_float_str.format(y_i),
                                                           "Longitude": ll_float_str.format(x_i),
                                                           "Depth": ll_float_str.format(z)})
        trace_element.append(loc_element)

    return trace_element


def fault_polygon_xml(polygon: Polygon, section_name: str, z: Union[float, int] = 0):
    """

    :param polygon: Should be lon lat
    :param section_name:
    :param z: Generally zero
    :return:
    """
    polygon_element = ElemTree.Element("ZonePolygon", attrib={"name": section_name})
    location_list = ElemTree.Element("LocationList")

    ll_float_str = "{:.4f}"
    x, y = polygon.exterior.xy
    for x_i, y_i in zip(x[:-1], y[:-1]):
        loc_element = ElemTree.Element("Location", attrib={"Latitude": ll_float_str.format(y_i),
                                                           "Longitude": ll_float_str.format(x_i),
                                                           "Depth": ll_float_str.format(z)})
        location_list.append(loc_element)
    polygon_element.append(location_list)
    return polygon_element
