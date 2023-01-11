"""
Classes that implement the Leapfrog fault model. Inherit from GenericFault and GenericMultiFault.
"""
from __future__ import annotations
import os
from typing import List, Union
from itertools import product, chain

import numpy as np
import geopandas as gpd
import pandas as pd

from shapely.affinity import translate
from shapely.ops import unary_union
from shapely.geometry import LineString, MultiLineString, Point, Polygon, MultiPoint

from fault_mesh.faults.generic import GenericMultiFault, GenericFault, normalize_bearing, smallest_difference
from fault_mesh.utilities.smoothing import smooth_trace
from fault_mesh.utilities.cutting import cut_line_between_two_points, cut_line_at_point
from fault_mesh.utilities.graph import connected_nodes, suggest_combined_name
from fault_mesh.faults.connected import ConnectedFaultSystem


class LeapfrogMultiFault(GenericMultiFault):
    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, sort_sr: bool = False,
                 segment_distance_tolerance: float = 100., smoothing_n_refinements: int = 5,
                 remove_colons: bool = True, tolerance: float = 100.):

        self._smoothing_n = smoothing_n_refinements
        super(LeapfrogMultiFault, self).__init__(fault_geodataframe=fault_geodataframe, 
                                                 sort_sr=sort_sr,
                                                 remove_colons=remove_colons)

        self._segment_distance_tolerance = segment_distance_tolerance
        self._cutting_hierarchy = []
        self._hierarchy_dict = None
        self._includes_connected = False
        self._curated_faults = None
        self._connections = None
        self._neighbour_connections = None
        self.connected_faults = None
        self._segment_dict = None
        self._multi_segment_dict = None
        self._inter_fault_connections = None

    def add_fault(self, series: pd.Series, depth_type: str = "D90", remove_colons: bool = False,
                  tolerance: float = 100.):
        cfmFault = LeapfrogFault.from_series(series, parent_multifault=self,
                                             remove_colons=remove_colons, tolerance=tolerance)
        cfmFault.smoothed_trace = smooth_trace(cfmFault.nztm_trace, n_refinements=self.smoothing_n)
        self.faults.append(cfmFault)

    @property
    def smoothing_n(self):
        return self._smoothing_n

    @property
    def segment_distance_tolerance(self):
        return self._segment_distance_tolerance

    @segment_distance_tolerance.setter
    def segment_distance_tolerance(self, value: float):
        self._segment_distance_tolerance = value

    @property
    def cutting_hierarchy(self):
        return self._cutting_hierarchy

    @cutting_hierarchy.setter
    def cutting_hierarchy(self, hierarchy: List[str]):
        for fault in self.curated_faults:
            if fault.name not in hierarchy:
                raise NameError(f"{fault.name} not included in cutting hierarchy")
        for name in hierarchy:
            if name not in self.names:
                print(f"Warning: unrecognised fault({name}) in cutting hierarchy")

        self._cutting_hierarchy = hierarchy

    def trim_cutting_hierarchy(self):
        trimmed_hierarchy = []
        for fault in self.cutting_hierarchy:
            if fault in self.names:
                trimmed_hierarchy.append(fault)
        self.cutting_hierarchy = trimmed_hierarchy

    def write_cutting_hierarchy(self, out_file: str):
        with open(out_file, "w") as out_id:
            for name in self.cutting_hierarchy:
                out_id.write(name + "\n")

    @property
    def hierarchy_dict(self):
        return self._hierarchy_dict

    @hierarchy_dict.setter
    def hierarchy_dict(self, dictionary):
        assert isinstance(dictionary, dict)
        self._hierarchy_dict = dictionary

    @property
    def curated_faults(self):
        if self._curated_faults is None:
            return self.faults
        else:
            return self._curated_faults

    @curated_faults.setter
    def curated_faults(self, fault_list: list):
        assert isinstance(fault_list, list)
        assert all([isinstance(fault, (ConnectedFaultSystem, LeapfrogFault)) for fault in fault_list])
        self._curated_faults = fault_list

    @property
    def curated_fault_dict(self):
        return {fault.name: fault for fault in self.curated_faults}

    def read_cutting_hierarchy(self, hierarchy_file: str):
        assert os.path.exists(hierarchy_file)
        with open(hierarchy_file, "r") as hfile:
            names = hfile.readlines()
            self.cutting_hierarchy = [name.strip() for name in list(names)]
        self.hierarchy_dict = {name: i for i, name in enumerate(self.cutting_hierarchy)}

    def suggest_cutting_hierarchy(self, prefix: str):
        name_ls = []
        sr_ls = []
        for fault in self.curated_faults:
            if isinstance(fault, LeapfrogFault):
                name_ls.append(fault.name)
                sr_ls.append(fault.sr_best)

            else:
                assert isinstance(fault, ConnectedFaultSystem)
                name_ls.append(fault.name)
                sr_ls.append(max([seg.sr_best for seg in fault.segments]))

        df = pd.DataFrame({"name": name_ls, "sr": sr_ls})
        df.sort_values(["sr", "name"], ascending=(False, True), inplace=True)

        out_names = list(df.name)
        out_file_name = prefix + "_suggested.csv"

        with open(out_file_name, "w") as out_id:
            for name in out_names:
                out_id.write(name + "\n")





    @property
    def names(self):
        return [fault.name for fault in self.curated_faults]

    def find_connections(self, verbose: bool = True):
        """
        Find all connections between faults in the fault list using networkx
        :param verbose: print out information about individual connections
        :return:
        """
        # Total number of connections of any type
        connections = []
        # Connections between faults that are neighbours (along strike), everything else is a branch
        neighbour_connections = []
        for fault in self.faults:
            for other_fault in self.faults:
                if other_fault.name != fault.name:
                    if fault.nztm_trace.distance(other_fault.nztm_trace) <= self.segment_distance_tolerance:
                        if verbose:
                            print(f"Connection: {fault.name} and {other_fault.name}")

                        connections.append([fault.name, other_fault.name])
                        # Check if faults are neighbours
                        conditions = []
                        for p1, p2 in product([fault.end1, fault.end2], [other_fault.end1, other_fault.end2]):
                            conditions.append(p1.distance(p2) <= self.segment_distance_tolerance)
                        if any(conditions):
                            neighbour_connections.append([fault.name, other_fault.name])


        self._connections = connections
        self._neighbour_connections = neighbour_connections
        print(f"Found {len(connections)} connections")
        print(f"Found {len(neighbour_connections)} connections between segment ends")

    @property
    def connections(self):
        if self._connections is None:
            self.find_connections()
        return self._connections

    @property
    def neighbour_connections(self):
        if self._neighbour_connections is None:
            self.find_connections()
        return self._neighbour_connections

    def suggest_fault_systems(self, out_prefix: str):
        connected = connected_nodes(self.neighbour_connections)
        out_name = out_prefix + "_suggested.csv"
        with open(out_name, "w") as out_id:
            for connected_set in connected:
                suggested_name = suggest_combined_name(connected_set)
                out_list = [suggested_name] + sorted(list(connected_set))
                out_str = ",".join(out_list) + "\n"
                out_id.write(out_str)

    def read_fault_systems(self, fault_system_csv: str):
        self.connected_faults = []
        with open(fault_system_csv, "r") as in_id:
            con_data = in_id.readlines()
            for line in con_data:
                elements = line.strip().split(",")
                name = elements[0].strip()
                name = name.replace(":", "")
                segs = [element.strip() for element in elements[1:]]
                if any([seg in self.names for seg in segs]):
                    cfault = ConnectedFaultSystem(overall_name=name, cfm_faults=self, segment_names=segs,
                                                  tolerance=self.segment_distance_tolerance)
                    self.connected_faults.append(cfault)

    def generate_curated_faults(self):
        curated_faults = []
        segments_already_included = []
        if len(self.connected_faults):
            curated_faults += self.connected_faults
            for fault in curated_faults:
                segments_already_included += [seg.name for seg in fault.segments]

        else:
            print("Warning: no multi-segment faults: is this what you wanted?")

        for fault in self.faults:
            if fault.name not in segments_already_included:
                curated_faults.append(fault)
                segments_already_included.append(fault.name)

        self.curated_faults = curated_faults

    @classmethod
    def from_nz_cfm_shp(cls, filename: str, exclude_region_polygons: List[Polygon] = None, depth_type: str = "D90",
                        exclude_region_min_sr: float = 1.8, include_names: list = None, exclude_aus: bool = True,
                        exclude_zero: bool = True, sort_sr: bool = False, remove_colons: bool = False,
                        smoothing_n_refinements: int = 5):

        trimmed_fault_gdf = cls.gdf_from_nz_cfm_shp(filename=filename, exclude_region_polygons=exclude_region_polygons,
                                                    depth_type=depth_type, exclude_region_min_sr=exclude_region_min_sr,
                                                    include_names=include_names, exclude_aus=exclude_aus,
                                                    exclude_zero=exclude_zero)
        multi_fault = cls(trimmed_fault_gdf, sort_sr=sort_sr,
                          remove_colons=remove_colons, smoothing_n_refinements=smoothing_n_refinements)

        return multi_fault

    @property
    def segment_dict(self):
        if self._segment_dict is None:
            self.generate_segment_dicts()
        return self._segment_dict

    @property
    def multi_segment_dict(self):
        if self._multi_segment_dict is None:
            self.generate_segment_dicts()
        return self._multi_segment_dict

    def generate_segment_dicts(self):
        segment_dict = {}
        multi_segment_dict = {}
        for fault in self.faults:
            segment_dict[fault.name] = fault
        self._segment_dict = segment_dict

        for fault in self.connected_faults:
            for seg in fault.segments:
                multi_segment_dict[seg.name] = fault
        self._multi_segment_dict = multi_segment_dict

    def find_inter_fault_connections(self):
        connection_set = set()
        for con in self.connections:
            f1 = self.segment_dict[con[0]]
            f2 = self.segment_dict[con[1]]

            f1_name = f1.name if not f1.is_segment else self.multi_segment_dict[f1.name].name
            f2_name = f2.name if not f2.is_segment else self.multi_segment_dict[f2.name].name
            if f1_name != f2_name:
                connection_set.add(tuple(sorted([f1_name, f2_name], key=lambda x: x.lower())))

        self._inter_fault_connections = connection_set

    @property
    def inter_fault_connections(self):
        if self._inter_fault_connections is None:
            self.find_inter_fault_connections()
        return self._inter_fault_connections

    def find_terminations(self, fault_name: str):
        fault_index = self.cutting_hierarchy.index(fault_name)
        terminations = []
        cons = [con for con in self.inter_fault_connections if fault_name in con]
        for con in cons:
            other_fault = [a for a in list(con) if a != fault_name][0]
            other_index = self.cutting_hierarchy.index(other_fault)
            if fault_index > other_index:
                terminations.append(con)

        return terminations


class LeapfrogFault(GenericFault):
    """
    Represents either a whole fault (for simple faults) or one segment. Behaviours is slightly
    """
    def __init__(self, parent_multifault: LeapfrogMultiFault = None, smoothing: int = 5,
                 trimming_gradient: float = 1.0, segment_distance_tolerance: float = 100.,
                 parent_connected=None):
        self._end1 = None
        self._end2 = None
        self._neighbouring_segments = None
        self._neighbour_dict = {}
        self._neighbour_angle_dict = {}

        super(LeapfrogFault, self).__init__(None)
        self._parent = parent_multifault

        self.connections = []
        self._is_segment = False
        self.parent_connected_fault = None
        self._smoothing = smoothing
        self._trimming_gradient = trimming_gradient
        self._segment_distance_tolerance = segment_distance_tolerance
        self._contours = None
        self._smoothed_trace = None
        self._footprint = None

    @property
    def is_segment(self):
        """
        Records whether instance is a segment of a larger multi-segment fault like the Alpine Fault.
        :return:
        """
        return self._is_segment

    @property
    def smoothing(self):
        """
        n value to use in Chaikin's corner-cutting algorithm.
        :return:
        """
        return self._smoothing

    @property
    def trimming_gradient(self):
        """
        Factor that controls how much the ends of segment contours of a multi-segment fault are
        shortened to allow
        :return:
        """

        return self._trimming_gradient

    @property
    def neighbouring_segments(self):
        return self._neighbouring_segments

    @neighbouring_segments.setter
    def neighbouring_segments(self, segment_list: list):
        assert isinstance(segment_list, list)
        if len(segment_list) > 2:
            raise ValueError(f"Too many ({len(segment_list):d}) neighbours supplied for segment {self.name:s}.\n"
                             f"Only two neighbours are allowed. Please turn extra neighbour(s) into separate faults.\n"
                             f"Neighbours supplied: {[neighbour.name for neighbour in segment_list]}")
        if len(segment_list) == 0:
            self._is_segment = False
        else:
            self._is_segment = True

        self._neighbouring_segments = segment_list

        for seg in segment_list:
            if seg.nztm_trace.distance(self.end1) <= self.segment_distance_tolerance:
                self._neighbour_dict[self.end1.coords[0]] = seg
                self._neighbour_angle_dict[self.end1.coords[0]] = self.neighbour_angle(seg)
            elif seg.nztm_trace.distance(self.end2) <= self.segment_distance_tolerance:
                self._neighbour_dict[self.end2.coords[0]] = seg
                self._neighbour_angle_dict[self.end2.coords[0]] = self.neighbour_angle(seg)
            else:
                raise ValueError(f"Orphan neighbour segment for {self.name}: {seg.name}")

    @property
    def neighbour_angle_dict(self):
        return self._neighbour_angle_dict

    @property
    def neighbour_dict(self):
        return self._neighbour_dict

    @property
    def strike(self):
        return normalize_bearing(self.dip_dir - 90.)

    @property
    def along_strike_vector(self):
        strike_rad = np.radians(self.strike)
        return np.array([np.sin(strike_rad), np.cos(strike_rad), 0.])

    @property
    def across_strike_vector(self):
        dip_dir_rad = np.radians(self.dip_dir)
        return np.array([np.sin(dip_dir_rad), np.cos(dip_dir_rad), 0.])

    @property
    def smoothed_trace(self):
        return self._smoothed_trace

    @smoothed_trace.setter
    def smoothed_trace(self, trace: LineString):
        assert isinstance(trace, LineString)
        self._smoothed_trace = trace

    @property
    def parent(self):
        """
        Return LeapfrogMultiFault instance that this fault is part of.
        :return:
        """
        return self._parent

    def depth_contour(self, depth: float, smoothing: bool = True, km=False):
        """
        Generate contour of fault surface at depth below surface
        :param depth: In metres, upwards is positive
        :param smoothing: N for use with Chaikin's corner cutting
        :param km: If True, divide depth by 1000
        :return: LineString or MultiLineString representing contour
        """
        if depth <= 0:
            shift = depth / self.down_dip_vector[-1]
        else:
            shift = (-1 * depth) / self.down_dip_vector[-1]
        if km:
            shift *= 1000.

        xo, yo, zo = shift * self.down_dip_vector

        if smoothing:
            assert self.smoothing is not None
            if self.is_segment:
                smoothed_contour = self.smoothed_trace

            else:
                self.smoothed_trace = smooth_trace(self.nztm_trace, self.smoothing)
                smoothed_contour = self.smoothed_trace


        else:
            smoothed_contour = self.nztm_trace

        shifted_contour = translate(smoothed_contour, xo, yo, zo)

        if self.is_segment:

            contour_e1 = Point(shifted_contour.coords[0])
            contour_e2 = Point(shifted_contour.coords[-1])

            flipping_conditions = [all([self.end1.x > self.end2.x, contour_e1.x < contour_e2.x]),
                                   all([self.end1.x < self.end2.x, contour_e1.x > contour_e2.x]),
                                   all([self.end1.x == self.end2.x, self.end1.y > self.end2.y,
                                        contour_e1.y < contour_e2.y]),
                                   all([self.end1.x == self.end2.x, self.end1.y < self.end2.y,
                                        contour_e1.y > contour_e2.y])]

            if any(flipping_conditions):
                contour_e2 = Point(shifted_contour.coords[0])
                contour_e1 = Point(shifted_contour.coords[-1])

            if len(self.neighbouring_segments) > 1:
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                e1_box_intersection = e1_box.boundary.intersection(self.nztm_trace)
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                e2_box_intersection = e2_box.boundary.intersection(self.nztm_trace)

                if not e1_box.intersects(e2_box):
                    trimmed_contour = cut_line_between_two_points(shifted_contour, [e1_box_intersection,
                                                                                     e2_box_intersection])
                else:
                    trimmed_contour = None

            elif self.end1.coords[0] in self.neighbour_dict.keys():
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                if all([point.within(e1_box) for point in (contour_e1, contour_e2)]):
                    trimmed_contour = None
                else:
                    e1_box_intersection = e1_box.boundary.intersection(self.nztm_trace)

                    if isinstance(e1_box_intersection, Point):
                        split_line = cut_line_at_point(shifted_contour, e1_box_intersection)

                        if split_line[0].distance(contour_e1) == 0.:
                            trimmed_contour = split_line[1]
                        elif split_line[1].distance(contour_e1) == 0.:
                            trimmed_contour = split_line[0]
                        else:
                            print("neither intersects")
                            trimmed_contour = None
                    else:
                        trimmed_contour = None

            else:  # self.end2 in self.neighbour_dict.keys()
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                if all([point.within(e2_box) for point in (contour_e1, contour_e2)]):
                    trimmed_contour = None
                else:
                    e2_box_intersection = e2_box.boundary.intersection(self.nztm_trace)
                    if isinstance(e2_box_intersection, Point):
                        split_line = cut_line_at_point(shifted_contour, e2_box_intersection)

                        if split_line[0].distance(contour_e2) == 0.:
                            trimmed_contour = split_line[1]
                        elif split_line[1].distance(contour_e2) == 0.:
                            trimmed_contour = split_line[0]
                        else:
                            print("neither intersects")
                            trimmed_contour = None
                    else:
                        trimmed_contour = None

            # if trimmed_contour is not None:
            #     if abs(trimmed_contour.length - smoothed_contour.length) < distance_tolerance:
            #         trimmed_contour = None

        else:
            trimmed_contour = shifted_contour

        return trimmed_contour

    def generate_depth_contours(self, depths: Union[np.ndarray, List[float]], smoothing: bool = True, damping: int = None,
                       km: bool = False):
        contours = [self.depth_contour(depth, smoothing, km) for depth in depths]

        if max(depths) > 0:
            depths *= -1

        self.contours =  gpd.GeoDataFrame({"depth": depths}, geometry=contours)

    @property
    def nztm_trace(self):
        return self._nztm_trace

    @nztm_trace.setter
    def nztm_trace(self, trace: LineString):
        assert isinstance(trace, (LineString, MultiLineString))
        if isinstance(trace, MultiLineString):
            trace = list(trace.geoms)[0]

        if trace.has_z:
            new_trace = LineString([(xi, yi, 0.) for xi, yi, _ in trace.coords])
        else:
            new_trace = LineString([(xi, yi, 0.) for xi, yi in trace.coords])
        self._nztm_trace = new_trace
        new_coord_array = np.array(new_trace.coords)
        self._end1 = Point(new_coord_array[0])
        self._end2 = Point(new_coord_array[-1])


    @property
    def end1(self):
        return self._end1

    @property
    def end2(self):
        return self._end2

    def clipping_box(self, centre_point: Point, along_half_width: float, across_half_width: float = 100000.):
        point_array = np.array(centre_point.coords)
        across_shift = across_half_width * self.across_strike_vector
        along_shift = along_half_width * self.along_strike_vector
        out_array = point_array + np.array([across_shift + along_shift,
                                            across_shift - along_shift,
                                            - 1 * (across_shift + along_shift),
                                            along_shift - across_shift])
        return Polygon(out_array)


    def end_clipping_box(self, end_i: Point, depth: float, gradient_adjustment: float = 1.,
                         across_half_width: float = 10000.):
        end_angle = self.neighbour_angle_dict[end_i.coords[0]]
        end_width = gradient_adjustment * np.tan(np.radians(end_angle)) * (depth / np.sin(np.radians(self.dip_best)))

        return self.clipping_box(end_i, end_width, across_half_width=across_half_width)

    def neighbour_angle(self, neighbour: LeapfrogFault):
        strike_diff = smallest_difference(self.dip_dir, neighbour.dip_dir)
        if strike_diff > 90.:
            strike_diff = 180. - strike_diff
            dip_diff = 180 - self.dip_best - neighbour.dip_best
            if dip_diff > 90.:
                raise ValueError(f"{self.name} and {neighbour.name}: shallow dips in opposite directions."
                                 "Are you sure you want to connect?")
        else:
            dip_diff = abs(neighbour.dip_best - self.dip_best)
        return max([2 * dip_diff, strike_diff]) / 2.

    @property
    def segment_distance_tolerance(self):
        if self.parent is not None:
            return self.parent.segment_distance_tolerance
        else:
            return self._segment_distance_tolerance

    @property
    def footprint(self):
        if self._footprint is None:
            self.calculate_footprint()
        return self._footprint


    @property
    def footprint_linestring(self):
        return LineString(self.footprint.exterior.coords)

    def calculate_footprint(self, smoothed: bool = True, buffer: float = 5000.):
        if smoothed:
            trace = self.smoothed_trace
        else:
            trace = self.nztm_trace

        buffer_top_offset = self.across_strike_vector * -1 * buffer
        shifted_top = translate(trace, *buffer_top_offset)

        bottom_trace = list(self.contours.geometry)[-1]
        buffer_bottom_offset = self.across_strike_vector * buffer
        shifted_bottom = translate(bottom_trace, *buffer_bottom_offset)

        buffer_combined = MultiLineString([shifted_top, shifted_bottom])
        self._footprint = buffer_combined.minimum_rotated_rectangle

    def find_terminations(self):
        return self.parent.find_terminations(self.name)

    def adjust_footprint(self):
        terms = list(set(chain(*self.find_terminations())))
        if terms:
            cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
            for fault in cutting_faults:
                for nearest_end, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
                    if nearest_end.distance(fault.nztm_trace) < 1.e3:
                        if isinstance(fault, ConnectedFaultSystem):
                            closest_seg = min(fault.segments, key=lambda x: x.nztm_trace.distance(nearest_end))
                        else:
                            closest_seg = fault

                        self.extend_footprint(nearest_end, other_end, closest_seg)

            # fp_to_merge = [self.footprint] + [fault.footprint for fault in cutting_faults]
            # merged_footprints = unary_union(fp_to_merge)
            # cutting_lines = []
            # for end_i, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
            #     if not any([end_i.distance(fault.nztm_trace) < self.segment_distance_tolerance for fault in cutting_faults]):
            #         cutting_lines.append((LineString([np.array(end_i) + self.across_strike_vector * line_length,
            #                                           np.array(end_i) - self.across_strike_vector * line_length]),
            #                               other_end))
            # if len(cutting_lines):
            #     if len(cutting_lines) > 1:
            #         print(f"{self.name}: more than one cutting line, choosing first...")
            #     splitter, other_end = cutting_lines[0]
            #     split_footprint = split(merged_footprints, splitter)
            #     kept_polys = [poly for poly in list(split_footprint) if other_end.within(poly)]
            #     if len(kept_polys) > 1:
            #         print(f"{self.name}: more than one cut polygon, choosing first...")
            #     self._footprint = kept_polys[0]
            #
            # else:
            #     self._footprint = merged_footprints

    def extend_footprint(self, end_i: Point, other_end: Point, other_segment: LeapfrogFault,
                         deepest_contour_depth: float = 30.e3, search_line_length: float = 1.5e5,
                         buffer_size: float = 5.e3, fall_back_distance: float = 40.3):
        """

        :param end_i: End to extend
        :param other_end: Other end of segment
        :param other_segment: Other
        :param deepest_contour_depth:
        :param search_line_length:
        :param buffer_size:
        :param fall_back_distance:
        :return:
        """
        # Find strike direction
        diff_vector = np.array(other_end.coords) - np.array(end_i.coords)
        if np.dot(diff_vector, self.along_strike_vector) > 0:
            strike_direction = -1 * self.along_strike_vector
        else:
            strike_direction = self.along_strike_vector


        if self.is_segment:
            bottom_trace = None
            contour_depth = deepest_contour_depth
            while bottom_trace is None:
                bottom_trace = self.depth_contour(contour_depth)
                contour_depth -= 2000.
        else:
            bottom_trace = list(self.contours.geometry)[-1]
        bot1 = np.array(bottom_trace.coords)[0]
        bot2 = np.array(bottom_trace.coords)[-1]

        def distance_along_strike(point: np.array, reference: np.array, strike_vector: np.array):
            diff = point - np.array(reference)
            return np.dot(diff, strike_vector)

        bot_i = bot1 if (distance_along_strike(bot1, np.array(end_i.coords), strike_direction) >
                         distance_along_strike(bot2, np.array(end_i.coords), strike_direction)) else bot2

        search_line = LineString([bot_i,
                                  bot_i + search_line_length * strike_direction])

        if other_segment.is_segment:
            other_fault = other_segment.parent_connected_fault
        else:
            other_fault = other_segment

        other_contour = other_fault.depth_contour(deepest_contour_depth)
        if other_contour.intersects(bottom_trace):
            return

        else:
            trace_intersection = search_line.intersection(other_fault.nztm_trace)
            if isinstance(trace_intersection, MultiPoint):
                trace_intersection = min(list(trace_intersection), key=lambda x: x.distance(Point(bot_i)))
            contour_intersection = search_line.intersection(other_contour)
            if isinstance(contour_intersection, MultiPoint):
                contour_intersection = min(list(contour_intersection), key=lambda x: x.distance(Point(bot_i)))

            if any([a.is_empty for a in [trace_intersection, contour_intersection]]):
                if all([a.is_empty for a in [trace_intersection, contour_intersection]]):
                    corner_point = Point(bot_i + fall_back_distance * strike_direction)
                elif not trace_intersection.is_empty:
                    corner_point = trace_intersection
                else:
                    corner_point = contour_intersection
            else:
                corner_point = max([contour_intersection, trace_intersection], key=lambda x: x.distance(Point(bot_i)))

            triangle = Polygon(np.vstack([bot_i, np.array(corner_point.coords), np.array(end_i.coords)])).buffer(buffer_size,
                                                                                        cap_style=2)


            if self.is_segment:
                new_boundary = unary_union([self.parent_connected_fault.footprint, triangle])
                self.parent_connected_fault._footprint = new_boundary
            else:
                new_boundary = unary_union([self.footprint, triangle])
                self._footprint = new_boundary

            return
