from typing import Union, List
import xml.etree.ElementTree as ElemTree
from xml.dom import minidom
import os

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import LineString, Polygon, MultiPolygon, MultiLineString
from shapely.ops import unary_union
from pyproj import Transformer
# import warnings

import logging

# from src.eq_fault_geom.geomio.cfm_logger import CfmLogger

transformer = Transformer.from_crs(2193, 4326, always_xy=True)
transform_wgs2nztm = Transformer.from_crs(4326, 2193, always_xy=True)

dip_direction_ranges = {"E": (45, 135), "NE": (0, 90), "N": (315, 45), "NW": (270, 360), "SW": (180, 270),
                        "S": (135, 225), "SE": (90, 180), "W": (225, 315)}
valid_dip_directions = list(dip_direction_ranges.keys()) + [None, "SUBVERTICAL AND VARIABLE"]

dominant_rake_ranges = {"reverse": (45, 135), "dextral": (135, 225), "sinistral": (315, 45), "normal": (225, 315)}
secondary_rake_ranges = {"reverse": (0, 180), "dextral": (90, 270), "sinistral": (270, 90), "normal": (180, 360)}
possible_rake_dirs = ['dextral', 'normal', 'reverse', 'sinistral', 'dextral and reverse', 'normal and dextral']

# List of subduction zones to exclude if desired
subduction_names = ("hikurangi", "puysegur")

dom_dict = {26: 'Puysegur subduction front',
            10: 'Hikurangi subduction front',
            20: 'Otago',
            6: 'East Cape block',
            2: 'North Waikato - South Auckland',
            21: 'Southern South Island',
            25: 'Puysegur Ridge - Bank',
            9: 'Hikurangi outer rise',
            14: 'Marlborough Fault System',
            12: 'Kapiti-Manawatu',
            7: 'Hikurangi accretionary margin',
            19: 'Southern Alps',
            23: 'Fiordland block',
            15: 'Alpine Fault',
            18: 'Central Canterbury',
            13: 'Northwestern South Island',
            4: 'Havre Trough - Taupo Rift',
            17: 'Northeastern Canterbury',
            24: 'Western Fiordland Margin - Caswell High',
            3: 'western North Island',
            16: 'North Mernoo fracture zone',
            27: 'Puysegur-Caswell High outer rise',
            5: 'North Island Dextral Fault Belt',
            8: 'Hikurangi accretionary margin - eastern zone'}

valid_dip_range = [0, 90]
valid_depth_range = [0, 50]
valid_rake_range = [0, 360]
valid_sr_range = [0, 60]

# These fields aren't crucial but are in some versions of the relevant files
expected_fields = ['D90', 'Depth_max', 'Depth_min', 'Dip_pref',
                   'Dip_dir', 'Dip_max', 'Dip_min', 'Name',
                   'Qual_Code', 'Rake_pref', 'Rake_plus', 'Rake_minus', 'Dom_sense',
                   'Sub_sense', 'Source1_1', 'Source2', 'SR_pref', 'SR_Max', 'SR_Min',
                   'geometry']

# There will be a mess if these fields don't exist
required_fields = ['Name', 'Fault_ID', 'geometry']


def decimal_deg_to_minutes(decdeg: float):
    """
    For use with hybrid model, which uses degrees and minutes.
    """
    if decdeg < 0.:
        decdeg *= -1
    deg = np.int(np.floor(decdeg))
    mins = 60. * (decdeg - deg)

    return deg, mins


def smallest_difference(value1, value2):
    """
    Finds smallest angle between two bearings
    :param value1:
    :param value2:
    :return:
    """
    abs_diff = abs(value1 - value2)
    if abs_diff > 180:
        smallest_diff = 360 - abs_diff
    else:
        smallest_diff = abs_diff

    return smallest_diff


def normalize_bearing(bearing: Union[float, int]):
    """
    change a bearing (in degrees) so that it is an azimuth between 0 and 360.
    :param bearing:
    :return:
    """
    while bearing < 0:
        bearing += 360.

    while bearing >= 360.:
        bearing -= 360.

    return bearing


def bearing_leq(value: Union[int, float], benchmark: Union[int, float], tolerance: Union[int, float] = 0.1):
    """
    Check whether a bearing (value) is anticlockwise of another bearing (benchmark)
    :param value:
    :param benchmark:
    :param tolerance: to account for rounding errors etc
    :return:
    """
    smallest_diff = smallest_difference(value, benchmark)
    if smallest_diff > tolerance:
        compare_value = normalize_bearing(value + smallest_diff)
        return abs(compare_value - normalize_bearing(benchmark)) <= tolerance
    else:
        return False


def bearing_geq(value: Union[int, float], benchmark: Union[int, float], tolerance: Union[int, float] = 0.1):
    """
    Check whether a bearing (value) is clockwise of another bearing (benchmark)
    :param value:
    :param benchmark:
    :param tolerance: to account for rounding errors etc
    :return:
    """
    smallest_diff = smallest_difference(value, benchmark)
    if smallest_diff > tolerance:
        compare_value = normalize_bearing(value - smallest_diff)
        return abs(compare_value - normalize_bearing(benchmark)) <= tolerance
    else:
        return False


def reverse_bearing(bearing: Union[int, float]):
    """
    180 degrees from supplied bearing
    :param bearing:
    :return:
    """
    assert isinstance(bearing, (float, int))
    assert 0. <= bearing <= 360.
    new_bearing = bearing + 180.

    # Ensure strike is between zero and 360 (bearing)
    return normalize_bearing(new_bearing)


def reverse_line(line: LineString):
    """
    Change the order that points in a LineString object are presented.
    Updated to work with 3d lines (has_z), September 2021
    Important for OpenSHA, I think
    :param line:
    :return:
    """
    assert isinstance(line, LineString)
    if line.has_z:
        x, y, z = np.array(line.coords).T
    else:
        x, y = np.array(line.coords).T
    x_back = x[-1::-1]
    y_back = y[-1::-1]

    if line.has_z:
        z_back = z[-1::-1]
        new_line = LineString([[xi, yi, zi] for xi, yi, zi in zip(x_back, y_back, z_back)])
    else:
        new_line = LineString([[xi, yi] for xi, yi in zip(x_back, y_back)])
    return new_line


def calculate_dip_direction(line: LineString):
    """
    Calculate the strike of a shapely linestring object with coordinates in NZTM,
    then adds 90 to get dip direction. Dip direction is always 90 clockwise from strike of line.
    :param line: Linestring object
    :return:
    """
    # Get coordinates
    x, y = line.xy
    x, y = np.array(x), np.array(y)
    # Calculate gradient of line in 2D
    px = np.polyfit(x, y, 1, full=True)
    gradient_x = px[0][0]

    if len(px[1]):
        res_x = px[1][0]
    else:
        res_x = 0

    py = np.polyfit(y, x, 1, full=True)
    gradient_y = py[0][0]
    if len(py[1]):
        res_y = py[1][0]
    else:
        res_y = 0

    if res_x <= res_y:
        # Gradient to bearing
        bearing = 180 - np.degrees(np.arctan2(gradient_x, 1))
    else:
        bearing = 180 - np.degrees(np.arctan2(1 / gradient_y, 1))
    bearing_vector = np.array([np.sin(np.radians(bearing)), np.cos(np.radians(bearing))])

    # Determine whether line object fits strike convention
    relative_x = x - x[0]
    relative_y = y - y[0]

    distances = np.matmul(np.vstack((relative_x, relative_y)).T, bearing_vector)
    num_pos = np.count_nonzero(distances >= 0)
    num_neg = np.count_nonzero(distances < 0)

    if num_neg > num_pos:
        bearing += 180.

    dip_direction = bearing
    # Ensure strike is between zero and 360 (bearing)
    while dip_direction < 0:
        dip_direction += 360.

    while dip_direction >= 360.:
        dip_direction -= 360.

    return dip_direction


def root_mean_square(value_array: Union[np.ndarray, list, tuple]):
    """
    Helper function to turn max and min to stdev for inclusion in XML.
    :param value_array: Differences of values (e.g. sr_min and sr_max) from mean.
    :return:
    """
    data_array = np.array(value_array)
    assert all([data_array.size > 0, data_array.ndim == 1])
    rms = np.sqrt(np.mean(np.square(data_array - np.mean(data_array))))
    return rms


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


class CfmMultiFault:
    """
    Class to hold data for multiple faults, read in from shapefile (and hopefully also tsurfaces)
    """

    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, exclude_region_polygons: list = None,
                 exclude_region_min_sr: float = 1.8, include_names: list = None, depth_type: str = "D90",
                 exclude_aus: bool = True, exclude_zero: bool = True, sort_sr: bool = False,
                 remove_colons: bool = False, min_dfc: float = 0.0):
        self.logger = logging.getLogger('cmf_logger')
        self.check_input1(fault_geodataframe)
        self.check_input2(fault_geodataframe)

        self._faults = []

        # If appropriate, clip out data that fall within exclude_regions
        if exclude_region_polygons is not None:
            assert isinstance(exclude_region_polygons, list)
            assert all([isinstance(a, Polygon) for a in exclude_region_polygons])
            exclude_regions_nztm = []
            # Check that polygons are in NZTM, otherwise convert them
            for poly in exclude_region_polygons:
                x, y = poly.exterior.xy
                if all(np.array(y) < 0):
                    # Assume in WGS (Lon Lat), convert to NZTM
                    new_x, new_y = transform_wgs2nztm.transform(np.array(x), np.array(y))
                    new_poly = Polygon([(xi, yi) for xi, yi in zip(new_x, new_y)])
                    exclude_regions_nztm.append(new_poly)
                else:
                    # Assume NZTM, do nothing
                    exclude_regions_nztm.append(poly)

            # Make list of faults outside region
            trimmed_fault_ls = []
            for i, row in fault_geodataframe.iterrows():
                if not any([row.geometry.within(poly) for poly in exclude_regions_nztm]):
                    trimmed_fault_ls.append(row)
                elif row["SR_pref"] >= exclude_region_min_sr:
                    trimmed_fault_ls.append(row)
                elif include_names is not None:
                    if row["Name"] in include_names:
                        trimmed_fault_ls.append(row)
                    else:
                        print(row["Name"])
                else:
                    print(row["Name"])
            trimmed_fault_gdf = gpd.GeoDataFrame(trimmed_fault_ls)

        else:
            trimmed_fault_gdf = fault_geodataframe

        # Temporarily avoid having to deal with zero dips
        trimmed_fault_gdf = trimmed_fault_gdf[trimmed_fault_gdf.Dip_pref > 0]

        # Exclude faults with zero slip rate
        if exclude_zero:
            trimmed_fault_gdf = trimmed_fault_gdf[trimmed_fault_gdf.SR_pref > 0.]

        # Exclude upper slope faults (A-US classification)
        if exclude_aus:
            trimmed_fault_gdf = trimmed_fault_gdf[trimmed_fault_gdf.Fault_stat != "A-US"]

        if depth_type is not "D90":
            assert min_dfc >= 0.
            trimmed_fault_gdf = trimmed_fault_gdf[trimmed_fault_gdf.Depth_Dfc >= min_dfc]

        if sort_sr:
            sorted_df = trimmed_fault_gdf.sort_values("SR_pref", ascending=False)
        else:
            # Sort alphabetically by name
            sorted_df = trimmed_fault_gdf.sort_values("Name")
        # Reset index to line up with alphabetical sorting
        sorted_df = sorted_df.reset_index(drop=True)
        for i, row in sorted_df.iterrows():
            self.add_fault(row, depth_type=depth_type, remove_colons=remove_colons)

        self.df = sorted_df

    def check_input1(self, fault_geodataframe):
        for field in required_fields:
            if field not in fault_geodataframe.columns:
                raise ValueError("Missing required field: {}".format(field))

    def check_input2(self, fault_geodataframe):
        for field in expected_fields:
            if field not in fault_geodataframe.columns:
                print("Warning: missing expected field: ({})".format(field))
                self.logger.warning("missing expected field")

    @property
    def faults(self):
        return self._faults

    def add_fault(self, series: pd.Series, depth_type: str = "D90", remove_colons: bool = False):
        cfmFault = CfmFault.from_series(series, parent_multifault=self, depth_type=depth_type,
                                        remove_colons=remove_colons)
        self.faults.append(cfmFault)

    @property
    def fault_numbers(self):
        if self.faults is not None:
            return [fault.number for fault in self.faults]
        else:
            return []

    @classmethod
    def from_shp(cls, filename: str, exclude_region_polygons: List[Polygon] = None, depth_type: str = "D90",
                 exclude_region_min_sr: float = 1.8, sort_sr: bool = False, min_dfc: float = 0.0):
        """
        Read CFM shapefile
        """
        assert os.path.exists(filename)
        fault_geodataframe = gpd.GeoDataFrame.from_file(filename)
        multi_fault = cls(fault_geodataframe, exclude_region_polygons=exclude_region_polygons,
                          exclude_region_min_sr=exclude_region_min_sr, depth_type=depth_type, sort_sr=sort_sr,
                          min_dfc=min_dfc)
        return multi_fault

    def to_opensha_xml(self, exclude_subduction: bool = True, buffer_width: float = 5000.,
                       write_buffers: bool = True):
        """
        Write out XML in OpenSHA format
        :param exclude_subduction: Do not include subduction zones from CFM
        :return:
        """
        assert self.faults
        assert isinstance(exclude_subduction, bool)
        # Base XML element
        opensha_element = ElemTree.Element("OpenSHA")
        # Fault model sub element
        fm_element = ElemTree.Element("FaultModel")
        opensha_element.append(fm_element)

        i = 0
        for fault in self.faults:
            # Identify subduction zone sources to include (only if exclude_subduction is True)
            exclude_condition = all([exclude_subduction, any([name in fault.name.lower()
                                                              for name in subduction_names])])
            # Add XML for fault
            if not exclude_condition:
                fm_element.append(fault.to_xml(section_id=i, buffer_width=buffer_width, write_buffers=write_buffers))
                i += 1

        # Awkward way of getting the xml file to be written in a way that's easy to read.
        # elmstr = ElemTree.tostring(opensha_element, encoding="UTF-8", xml_declaration=True)
        elmstr = ElemTree.tostring(opensha_element, encoding="UTF-8")
        xml_dom = minidom.parseString(elmstr)
        pretty_xml_str = xml_dom.toprettyxml(indent="  ", encoding="utf-8")

        return pretty_xml_str

    def to_xml_file(self, filename: str, exclude_subduction: bool = True):
        """
        Write to file
        :param filename:
        :param exclude_subduction: Do not include subduction zones from CFM
        :return:
        """

        with open(filename, "wb") as f:
            opnxml = self.to_opensha_xml(exclude_subduction=exclude_subduction)
            f.write(opnxml)

    def to_hybrid_csv(self, output_csv: str, exclude_subduction: bool = True):
        pass
        """
        :param output_csv: File to write output
        :param exclude_subduction: Do not include subduction zones from CFM
        """
        columns = ["FaultName", "DipMin", "DipBest", "DipMax", "DipDir", "DepthMin", "DepthBest", "DepthMax",
                   "Top", "SRMin", "SRBest", "SRMax"]

        output_str = ",".join(columns) + "\n"
        num_columns = len(columns)

        for fault in self.faults:
            # Identify subduction zone sources to include (only if exclude_subduction is True)
            exclude_condition = all([exclude_subduction, any([name in fault.name.lower()
                                                              for name in subduction_names])])
            # Add XML for fault
            if not exclude_condition:
                output_str += fault.fault_to_hybrid_csv(num_columns=num_columns)

        with open(output_csv, "w") as outfile:
            outfile.write(output_str)

    def to_gmt(self, output_csv: str, exclude_subduction: bool = True):
        """
        :param output_csv: File to write output
        :param exclude_subduction: Do not include subduction zones from CFM
        """
        out_str = ""
        for fault in self.faults:
            # Identify subduction zone sources to include (only if exclude_subduction is True)
            exclude_condition = all([exclude_subduction, any([name in fault.name.lower()
                                                              for name in subduction_names])])
            # Add XML for fault
            if not exclude_condition:
                out_str += fault.trace_to_gmt()

        with open(output_csv, "w") as out_file:
            out_file.write(out_str)


class CfmFault:
    def __init__(self, parent_multifault: CfmMultiFault = None):
        """

        :param parent_multifault:
        """
        # Attributes usually provided in CFM trace shapefile
        self.logger = logging.getLogger('cmf_logger')
        self._parent = parent_multifault
        self._depth_best, self._depth_min, self._depth_max, self._depth_stdev = (None,) * 4
        self._dip_best, self._dip_min, self._dip_max, self._dip_dir, self._dip_dir_str = (None,) * 5
        self._name = None
        self._number, self._qual_code = (None,) * 2
        self._rake_best, self._rake_max, self._rake_min = (None,) * 3
        self._sense_dom, self._sense_sec = (None,) * 2
        self._source1_1, self.source2 = (None,) * 2
        self._sr_best, self._sr_max, self._sr_min = (None,) * 3
        self._nztm_trace = None
        self._dom_num = None
        self._dom_name = None

        # Attributes required for OpenSHA XML
        self._section_id, self._section_name = (None,) * 2
        self._dip_sigma, self._rake_sigma, self._sr_sigma = (None,) * 3

    # Depths
    @property
    def depth_best(self):
        return self._depth_best

    @property
    def depth_max(self):
        return self._depth_max

    @property
    def depth_min(self):
        return self._depth_min

    @depth_best.setter
    def depth_best(self, depth: Union[float, int]):
        depth_v = self.validate_depth(depth)
        if self.depth_min is not None:
            if depth_v < self.depth_min:
                # print("{}: depth_best ({:.2f}) lower than depth_min ({:.2f})".format(self.name, depth_v,
                #                                                                       self.depth_min))
                self.logger.warning("depth_best lower than depth_min")

        if self.depth_max is not None:
            if depth_v > self.depth_max:
                # print("{}: depth_best ({:.2f}) greater than depth_max ({:.2f})".format(self.name, depth_v,
                #                                                                        self.depth_max))

                self.logger.warning("depth_best greater than depth_max ")
        self._depth_best = depth_v

    @depth_max.setter
    def depth_max(self, depth: Union[float, int]):
        depth_v = self.validate_depth(depth)
        for depth_value in (self.depth_min, self.depth_best):
            if depth_value is not None and depth_v < depth_value:
                print("Warning: depth_max lower than either depth_min or depth_best ({})".format(self.name))
                self.logger.warning("depth_max lower than either depth_min or depth_best")

        self._depth_max = depth_v

    @depth_min.setter
    def depth_min(self, depth: Union[float, int]):
        depth_v = self.validate_depth(depth)
        for depth_value in (self.depth_max, self.depth_best):
            if depth_value is not None and depth_v > depth_value:
                print("Warning: depth_min higher than either depth_max or depth_best ({})".format(self.name))
                self.logger.warning("depth_min higher than either depth_max or depth_best")
        self._depth_min = depth_v

    def validate_depth(self, depth: Union[float, int]):
        assert isinstance(depth, (float, int))
        depth_positive = depth if depth >= 0 else depth * -1
        if not valid_depth_range[0] < depth_positive <= valid_depth_range[1]:
            # print("{}: Supplied (lower) depth should be > {:.1f} and <= {:.1f}".format(self.name, valid_depth_range[0],
            #                                                                            valid_depth_range[1]))
            pass
        return depth_positive

    @property
    def depth_stdev(self):
        return self._depth_stdev

    @depth_stdev.setter
    def depth_stdev(self, dstd):
        assert isinstance(dstd, (int, float))
        assert dstd >= 0.
        self._depth_stdev = dstd

    # Dips
    @property
    def dip_best(self):
        return self._dip_best

    @property
    def dip_max(self):
        return self._dip_max

    @property
    def dip_min(self):
        return self._dip_min

    @property
    def dip_dir_str(self):
        return self._dip_dir_str

    @dip_best.setter
    def dip_best(self, dip: Union[float, int]):
        dip_v = self.validate_dip(dip)
        if self.dip_min is not None:
            if dip_v < self.dip_min:
                print("{}: dip_best ({}) lower than dip_min ({})".format(self.name, dip_v, self.dip_min))
        if self.dip_max is not None:
            if dip_v > self.dip_max:
                print("{}: dip_best ({}) greater than dip_max ({})".format(self.name, dip_v, self.dip_max))
        self._dip_best = dip_v

    @dip_max.setter
    def dip_max(self, dip: Union[float, int]):
        dip_v = self.validate_dip(dip)
        for key, dip_value in zip(["dip_min", "dip_best"], [self.dip_min, self.dip_best]):
            if dip_value is not None and bearing_leq(dip_v, dip_value):
                print("{}: dip_max ({}) is lower than {} ({})".format(self.name, dip_v, key, dip_value))
                self.logger.warning("dip_max is lower than dip min or dip best")
        self._dip_max = dip_v

    @dip_min.setter
    def dip_min(self, dip: Union[float, int]):
        dip_v = self.validate_dip(dip)
        for key, dip_value in zip(["dip_max", "dip_best"], [self.dip_max, self.dip_best]):
            if dip_value is not None and bearing_geq(dip_v, dip_value):
                # print("{}: dip_min ({:.2f}) is higher than {} ({:.2f})".format(self.name, dip_v, key, dip_value)) #{:.2f} formatting is not working
                print("{}: dip_min ({}) is higher than {} ({})".format(self.name, dip_v, key, dip_value))
                self.logger.warning("dip_min is higher than dip max or dip best")

        self._dip_min = dip_v

    @property
    def dip_sigma(self):
        if self._dip_sigma is not None:
            return self._dip_sigma
        elif not any([a is None for a in (self.dip_min, self.dip_max)]):
            return root_mean_square(np.array([self.dip_min, self.dip_max]))
        else:
            raise ValueError(
                "Insufficient data to calculate dip_sigma!")  # <= not sure if you need this step as this will go through validate_dip

    @dip_dir_str.setter
    def dip_dir_str(self, dip_dir: str):
        assert any([isinstance(dip_dir, str), dip_dir is None])
        if isinstance(dip_dir, str):
            if not dip_dir.upper() in valid_dip_directions:
                print("Unrecognised dip direction: {}".format(self.name))
                print("{}".format(dip_dir))
                self.validate_dip_direction()
            else:
                self._dip_dir_str = dip_dir.upper()
                if self.nztm_trace is not None:
                    self.validate_dip_direction()

        else:
            self._dip_dir_str = None
            if self.nztm_trace is not None:
                dd_from_trace = calculate_dip_direction(self.nztm_trace)
                self._dip_dir = dd_from_trace

    @property
    def dip_dir(self):
        """
        Azimuth (in degrees) of dip direction. Usually calculated from fault trace.
        :return:
        """
        return self._dip_dir

    def validate_dip_direction(self, tolerance: float = 10.):
        """
        Compares dip direction string (e.g. NW) with
        :return:
        """
        if any([a is None for a in [self.dip_dir_str, self.nztm_trace]]):
            print("Insufficient information to validate dip direction")
            if self.nztm_trace is not None and self._dip_dir is None:
                dip_dir = calculate_dip_direction(self.nztm_trace)
                self._dip_dir = dip_dir
            self.logger.warning("Insufficient information to validate dip direction")
            return
        else:
            # Trace and dip direction

            dd_from_trace = calculate_dip_direction(self.nztm_trace)
            if self.dip_dir_str != "SUBVERTICAL AND VARIABLE":
                min_dd_range, max_dd_range = dip_direction_ranges[self.dip_dir_str]
                if self.dip_dir_str != "N":
                    if not all([min_dd_range - tolerance <= dd_from_trace, dd_from_trace <= max_dd_range + tolerance]):
                        reversed_dd = reverse_bearing(dd_from_trace)
                        if all([min_dd_range - tolerance <= reversed_dd, reversed_dd <= max_dd_range + tolerance]):
                            self._nztm_trace = reverse_line(self.nztm_trace)
                            self._dip_dir = reversed_dd
                        else:
                            print("{}: Supplied trace and dip direction {} are inconsistent: expect either {:.1f}"
                                  "or {:.1f} dip azimuth. Please check...".format(self.name, self.dip_dir_str,
                                                                                  dd_from_trace, reversed_dd))
                            self.logger.warning("Supplied trace and dip direction are inconsistent")

                            self._dip_dir = dd_from_trace

                    else:
                        self._dip_dir = dd_from_trace
                else:
                    reversed_dd = reverse_bearing(dd_from_trace)
                    if any([315. - tolerance <= dds for dds in [dd_from_trace, reversed_dd]]):
                        self._dip_dir = max([dd_from_trace, reversed_dd])
                    elif any([dds <= 45. + tolerance for dds in [dd_from_trace, reversed_dd]]):
                        self._dip_dir = min([dd_from_trace, reversed_dd])
                    else:
                        print("{}: Supplied trace and dip direction {} are inconsistent: expect either {:.1f} or {:.1f}"
                              " dip azimuth. Please check...".format(self.name, self.dip_dir_str,
                                                                     dd_from_trace, reversed_dd))
                        self.logger.warning("Supplied trace and dip direction are inconsistent")

                        self._dip_dir = dd_from_trace

            else:
                self._dip_dir = dd_from_trace
            return

    @staticmethod
    def validate_dip(dip: Union[float, int]):
        """
        Generally between 0 and 90
        """
        assert isinstance(dip, (float, int))
        assert valid_dip_range[0] <= dip <= valid_dip_range[1]
        return dip

    @property
    def down_dip_vector(self):
        """
        Calculated from dip and dip direction
        """
        assert self.dip_best is not None
        if self.dip_dir is None:
            # Assume vertical
            return np.array([0., 0., -1])
        else:
            z = np.sin(np.radians(self.dip_best))
            x, y = np.cos(np.radians(self.dip_best)) * np.array([np.sin(np.radians(self.dip_dir)),
                                                                 np.cos(np.radians(self.dip_dir))])
        return np.array([x, y, -z])

    @property
    def down_dip_polygon(self):
        """

        """
        surface_trace_array = np.vstack(self.nztm_trace.xy).T
        if abs(self.down_dip_vector[-1]) > 1.e-3:
            bottom_trace = surface_trace_array + -1. * self.down_dip_vector[:-1] * \
                           (self.depth_best / self.down_dip_vector[-1]) * 1.e3
        else:
            bottom_trace = surface_trace_array + self.down_dip_vector[:-1] * 100.e3

        combined_polygon_array = np.vstack((surface_trace_array, bottom_trace[::-1]))
        return Polygon(combined_polygon_array)

    def surface_trace_buffer(self, buffer_distance: Union[int, float] = 100., cap_style: int = 2,
                             join_style: int = 2):
        """

        :param buffer_distance:
        :param cap_style: Default is flat,
        see https://shapely.readthedocs.io/en/stable/manual.html#shapely.geometry.CAP_STYLE
        :param join_style: Default is mitred
        :return:
        """
        return self.nztm_trace.buffer(buffer_distance, cap_style=cap_style, join_style=join_style)

    def combined_buffer_polygon(self, buffer_distance, cap_style: int = 2, join_style: int = 2, wgs: bool = True):
        trace_buffer = self.surface_trace_buffer(buffer_distance, cap_style=cap_style, join_style=join_style)
        combined_buffer = unary_union([trace_buffer, self.down_dip_polygon])

        if wgs:
            x, y = combined_buffer.exterior.xy
            wgs_x, wgs_y = transformer.transform(x, y)
            return Polygon([[xi, yi] for xi, yi in zip(wgs_x, wgs_y)])
        else:
            return combined_buffer

    # Trace
    @property
    def nztm_trace(self):
        return self._nztm_trace

    @nztm_trace.setter
    def nztm_trace(self, trace: LineString):
        if isinstance(trace, LineString):
            self._nztm_trace = trace
        else:
            assert isinstance(trace, MultiLineString)
            self._nztm_trace = list(trace.geoms)[0]

    @property
    def wgs_trace(self):
        if self.nztm_trace is not None:
            nztm_x, nztm_y = [np.array([a]).flatten() for a in self.nztm_trace.xy]
            wgs_x, wgs_y = transformer.transform(nztm_x, nztm_y)
            return LineString([(xi, yi) for xi, yi in zip(wgs_x, wgs_y)])

        else:
            return None

    # Rake and sense of slip
    @property
    def rake_best(self):
        return self._rake_best

    @property
    def rake_max(self):
        return self._rake_max

    @property
    def rake_min(self):
        return self._rake_min

    @property
    def sense_dom(self):
        return self._sense_dom

    @property
    def sense_sec(self):
        return self._sense_sec

    @rake_best.setter
    def rake_best(self, rake: Union[float, int]):
        rake_v = self.validate_rake(rake)
        if self.rake_min is not None:
            if rake_v < self.rake_min:
                # print("{}: rake_best ({.2f}) lower than rake_min ({.2f})".format(self.name, rake_v, self.rake_min))
                print("{}: rake_best ({}) lower than rake_min ({})".format(self.name, rake_v, self.rake_min))
                self.logger.warning("rake_best is lower than rake_min")
        if self.rake_max is not None:
            if rake_v > self.rake_max:
                # print("{}: rake_best ({.2f}) greater than rake_max ({.2f})".format(self.name, rake_v, self.rake_max))
                print("{}: rake_best ({}) greater than rake_max ({})".format(self.name, rake_v, self.rake_max))
                self.logger.warning("rake_best is greater than rake_max")
        self._rake_best = rake_v

        if self.sense_dom is not None:
            self.validate_rake_sense()

    @staticmethod
    def rake_to_opensha(rake: Union[float, int]):
        """
        To give opensha convention
        """
        new_rake = rake
        while new_rake > 180:
            new_rake -= 360.
        while new_rake <= -180.:
            new_rake += 360.
        return new_rake

    @rake_max.setter
    def rake_max(self, rake: Union[float, int]):
        rake_v = self.validate_rake(rake)
        for key, rake_value in zip(["rake_min", "rake_best"], [self.rake_min, self.rake_best]):
            if rake_value is not None and bearing_leq(rake_v, rake_value):
                print("{}: rake_max ({}) is lower than {} ({})".format(self.name, rake_v, key, rake_value))
                self.logger.warning("rake_max is lower than rake min or rake best")
        self._rake_max = rake_v

    @rake_min.setter
    def rake_min(self, rake: Union[float, int]):
        rake_v = self.validate_rake(rake)
        for key, rake_value in zip(["rake_max", "rake_best"], [self.rake_max, self.rake_best]):
            if rake_value is not None and bearing_geq(rake_v, rake_value):
                # print("{}: rake_min ({:.2f}) is higher than {} ({:.2f})".format(self.name, rake_v, key, rake_value))
                print("{}: rake_min ({}) is higher than {} ({})".format(self.name, rake_v, key, rake_value))
                self.logger.warning("rake_min is higher than rake max or rake best")
        self._rake_min = rake_v

    @staticmethod
    def validate_rake(rake: Union[float, int]):
        assert isinstance(rake, (float, int))
        if -180. <= rake < 0.:
            rake += 360
        while rake >= 360.:
            rake -= 360.
        assert valid_rake_range[0] <= rake <= valid_rake_range[1]
        return rake

    @sense_dom.setter
    def sense_dom(self, sense: str):
        assert any([isinstance(sense, str), sense is None])
        if sense is None:
            print("{}: Unexpected sense_dom: {}".format(self.name, sense))
        elif sense.lower() not in possible_rake_dirs:
            print("{}: Unexpected sense_dom: {}".format(self.name, sense.lower()))
        self._sense_dom = sense

    @sense_sec.setter
    def sense_sec(self, sense: str):
        assert any([sense is None, isinstance(sense, str)])
        if sense is not None:
            if sense.lower() not in possible_rake_dirs:
                print("{}: Unexpected sense_sec: {}".format(self.name, sense.lower()))
            if self.rake_best is not None:
                self.validate_rake_sense()
        self._sense_sec = sense

    def validate_rake_sense(self):
        if any([a is None for a in (self.rake_best, self.sense_dom)]):
            print("{}: Insufficient data to compare rake and slip sense".format(self.name))
            self.logger.warning("Insufficient data to compare rake and slip sense")
            return
        else:
            dominant_range = dominant_rake_ranges[self.sense_dom]
            if not all([bearing_geq(self.rake_best, dominant_range[0]),
                        bearing_leq(self.rake_best, dominant_range[1])]):
                print("{}: Supplied rake ({:.2f} deg) differs from dominant slip sense ({})".format(self.name,
                                                                                                    self.rake_best,
                                                                                                    self.sense_dom))
                self.logger.warning("Supplied rake differs from dominant slip sense")
            if self.sense_sec is not None:
                sec_range = secondary_rake_ranges[self.sense_sec]
                if not all([bearing_geq(self.rake_best, sec_range[0]),
                            bearing_leq(self.rake_best, sec_range[1])]):
                    print("{}: Supplied rake ({:.2f} deg) inconsistent with sec slip sense ({})".format(self.name,
                                                                                                        self.rake_best,
                                                                                                        self.sense_sec))
                    self.logger.warning("Supplied rake inconsistent with sec slip sense")

    @property
    def sr_best(self):
        return self._sr_best

    @property
    def sr_min(self):
        return self._sr_min

    @property
    def sr_max(self):
        return self._sr_max

    @sr_best.setter
    def sr_best(self, slip_rate: Union[float, int]):
        slip_rate = self.validate_sr(slip_rate)
        if self.sr_min is not None:
            if slip_rate < self.sr_min:
                print("{}: sr_best ({.2f}) lower than sr_min ({.2f})".format(self.name, slip_rate, self.sr_min))
        if self.sr_max is not None:
            if slip_rate > self.sr_max:
                print("{}: sr_best ({.2f}) greater than sr_max ({.2f})".format(self.name, slip_rate, self.sr_max))
        self._sr_best = slip_rate

    @sr_min.setter
    def sr_min(self, slip_rate: Union[float, int]):
        slip_rate = self.validate_sr(slip_rate)
        for key, sr_value in zip(["sr_max", "sr_best"], [self.sr_max, self.sr_best]):
            if sr_value is not None and bearing_geq(slip_rate, sr_value):
                print("{}: sr_min ({:.2f}) is higher than {} ({:.2f})".format(self.name, slip_rate, key, sr_value))
        self._sr_min = slip_rate

    @sr_max.setter
    def sr_max(self, slip_rate: Union[float, int]):
        self.validate_sr(slip_rate)
        for key, sr_value in zip(["sr_min", "sr_best"], [self.sr_max, self.sr_best]):
            if sr_value is not None and bearing_leq(slip_rate, sr_value):
                print("{}: sr_max ({:.2f}) is lower than {} ({:.2f})".format(self.name, slip_rate, key, sr_value))
        self._sr_max = slip_rate

    @staticmethod
    def validate_sr(slip_rate: Union[float, int]):
        assert isinstance(slip_rate, (float, int))
        assert valid_sr_range[0] <= slip_rate <= valid_sr_range[1]
        return slip_rate

    @property
    def sr_sigma(self):
        if self._sr_sigma is not None:
            return self._sr_sigma
        elif not any([a is None for a in (self.sr_min, self.sr_max)]):
            return root_mean_square(np.array([self.sr_min, self.sr_max]))
        else:
            print("{}: Insufficient data to calculate sr_sigma!".format(self.name))
            print(self.sr_min, self.sr_best, self.sr_max)
            return 0

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name_str: str):
        assert any([isinstance(name_str, str), name_str is None])
        if not name_str:
            print("Warning: empty name")
            name_str = "None"
        self._name = name_str

    @property
    def number(self):
        return self._number

    @number.setter
    def number(self, fault_number):
        assert isinstance(fault_number, int)
        if self.parent is not None:
            if fault_number in self.parent.fault_numbers:
                print("Duplicate fault number: {:d}".format(fault_number))

    @property
    def dom_name(self):
        return self._dom_name

    @dom_name.setter
    def dom_name(self, name):
        assert name in dom_dict.values()
        if self.dom_num is not None:
            assert dom_dict[self.dom_num] == name

        self._dom_name = name

    @property
    def dom_num(self):
        return self._dom_num

    @dom_num.setter
    def dom_num(self, num: int):
        assert isinstance(num, int)
        assert num in dom_dict.keys()
        self._dom_num = num

    @property
    def parent(self):
        return self._parent

    @classmethod
    def from_series(cls, series: pd.Series, parent_multifault: CfmMultiFault = None, depth_type: str = "D90",
                    remove_colons: bool = False):
        assert isinstance(series, pd.Series)
        assert depth_type in ["D90", "Dfc80", "Dfc66", "Dfc"]
        fault = cls(parent_multifault=parent_multifault)
        fault.name = series["Name"]
        if remove_colons:
            fault.name = series["Name"].replace(":", "")
        fault.number = int(series["Fault_ID"])
        fault.dip_best, fault.dip_min, fault.dip_max = series["Dip_pref"], series["Dip_min"], series["Dip_max"]
        fault.nztm_trace = series["geometry"]
        fault.dip_dir_str = series["Dip_dir"]
        fault.rake_best, fault.rake_min, fault.rake_max = series["Rake_pref"], series["Rake_minus"], series["Rake_plus"]
        fault.sense_dom, fault.sense_sec = series["Dom_sense"], series["Sub_sense"]
        fault.sr_best, fault.sr_min, fault.sr_max = series["SR_pref"], series["SR_min"], series["SR_max"]
        if depth_type == "D90":
            fault.depth_best = series["Depth_D90"]
        elif depth_type == "Dfc80":
            fault.depth_best = series["Depth_Dfc"] * 0.8
        elif depth_type == "Dfc66":
            fault.depth_best = series["Depth_Dfc"] * 0.666
        else:
            fault.depth_best = series["Depth_Dfc"]

        fault.dom_num = series["Domain_No"]
        fault.dom_name = series["DomainName"]

        return fault

    @classmethod
    def from_series_old(cls, series: pd.Series, parent_multifault: CfmMultiFault = None):
        assert isinstance(series, pd.Series)
        fault = cls(parent_multifault=parent_multifault)
        fault.name = series["FZ_Name"]
        fault.number = int(series["Number"])
        fault.dip_best, fault.dip_min, fault.dip_max = series["Dip_Best"], series["Dip_Min"], series["Dip_Max"]
        fault.nztm_trace = series["geometry"]
        fault.dip_dir_str = series["Dip_Dir"]
        fault.rake_best, fault.rake_min, fault.rake_max = series["Rake_Best"], series["Rake_Min"], series["Rake_Max"]
        fault.sense_dom, fault.sense_sec = series["Sense_Dom"], series["Sense_Sec"]
        fault.sr_best, fault.sr_min, fault.sr_max = series["SR_Best"], series["SR_Min"], series["SR_Max"]
        fault.depth_best, fault.depth_min, fault.depth_max = series["Depth_Best"], series["Depth_Min"], series[
            "Depth_Max"]
        return fault

    def to_xml(self, section_id: int, buffer_width: float = 5000., write_buffers: bool = True):
        # Unique fault identifier
        tag_name = "i{:d}".format(section_id)
        # Metadata
        attribute_dic = {"sectionId": "{:d}".format(section_id),
                         "sectionName": self.name,
                         "aveLongTermSlipRate": "{:.2f}".format(self.sr_best),
                         "slipRateStdDev": "{:.2f}".format(self.sr_sigma),
                         "aveDip": "{:.1f}".format(self.dip_best),
                         "aveRake": "{:.1f}".format(self.rake_to_opensha(self.rake_best)),
                         "aveUpperDepth": "0.0",
                         "aveLowerDepth": "{:.1f}".format(self.depth_best),
                         "aseismicSlipFactor": "0.0",
                         "couplingCoeff": "1.0",
                         "dipDirection": "{:.1f}".format(self.dip_dir),
                         "parentSectionId": "-1",
                         "connector": "false",
                         "domainNo": f"{self.dom_num:d}",
                         "domainName": f"{self.dom_name}"
                         }

        # Initialize XML element
        fault_element = ElemTree.Element(tag_name, attrib=attribute_dic)
        # Add sub element for fault trace
        trace_element = fault_trace_xml(self.wgs_trace, self.name)
        fault_element.append(trace_element)
        if write_buffers:
            # Add sub element for FZ buffer
            polygon_element = fault_polygon_xml(self.combined_buffer_polygon(buffer_width), self.name)
            fault_element.append(polygon_element)

        return fault_element

    def trace_to_gmt(self):
        x, y = self.wgs_trace.xy
        out_str = f">{self.name}\n"
        for xi, yi in zip(x, y):
            out_str += f"{xi:.6f} {yi:.6f}\n"
        return out_str

    @staticmethod
    def pad_commas(in_string: str, num_columns: int):
        num_commas = in_string.count(",")
        out_str = in_string.strip() + (num_columns - num_commas - 1) * "," + "\n"
        return out_str

    def to_segment(self, x0: float, y0: float, x1: float, y1: float, num_columns: int):
        x0_deg, x0_mins = decimal_deg_to_minutes(x0)
        y0_deg, y0_mins = decimal_deg_to_minutes(y0)
        x1_deg, x1_mins = decimal_deg_to_minutes(x1)
        y1_deg, y1_mins = decimal_deg_to_minutes(y1)

        seg_str = (f"{self.name},{x0_deg:d},{x0_mins:.1f},{y0_deg:d},{y0_mins:.1f},"
                   f"{x1_deg:d},{x1_mins:.1f},{y1_deg:d},{y1_mins:.1f}\n")
        return self.pad_commas(seg_str, num_columns)

    def trace_to_hybrid_csv(self, num_columns: int):
        """

        """

        x, y = self.wgs_trace.xy
        coord_ls = list(self.wgs_trace.coords)

        out_str = self.pad_commas(f"{self.name},{len(coord_ls) - 1}D", num_columns=num_columns)
        out_str += self.to_segment(x[0], y[0], x[-1], y[-1], num_columns=num_columns)

        for c0, c1 in zip(coord_ls[:-1], coord_ls[1:]):
            out_str += self.to_segment(c0[0], c0[1], c1[0], c1[1], num_columns=num_columns)

        out_str += self.pad_commas(f"{self.name},-1\n", num_columns=num_columns)
        return out_str

    def fault_to_hybrid_csv(self, num_columns: int):
        data_str = (f"{self.name},{int(self.dip_min):d},{int(self.dip_best):d},{int(self.dip_max):d},"
                    f"{self.dip_dir:.1f},{self.depth_best - self.depth_stdev:.1f},{self.depth_best:.1f},"
                    f"{self.depth_best + self.depth_stdev:.1f},0.0,{self.sr_min:.2f},{self.sr_best:.2f},"
                    f"{self.sr_max:.2f}\n")
        data_str += self.trace_to_hybrid_csv(num_columns=num_columns)
        return data_str
