from __future__ import annotations
import os
from typing import List, Union
from itertools import product, chain

import numpy as np
from eq_fault_geom.geomio.cfm_faults import CfmMultiFault, CfmFault, normalize_bearing, smallest_difference
import geopandas as gpd
import pandas as pd

from shapely.affinity import translate
from shapely.ops import unary_union, split
from shapely.geometry import LineString, MultiLineString, Point, Polygon

from fault_mesh.smoothing import straighten, smooth_trace
from fault_mesh.utilities.cutting import cut_line_between_two_points, cut_line_at_point
from fault_mesh.utilities.graph import connected_nodes, suggest_combined_name
from fault_mesh.connections import ConnectedFaultSystem


class LeapfrogMultiFault(CfmMultiFault):
    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, exclude_region_polygons: list = None,
                 exclude_region_min_sr: float = 1.8, include_names: list = None, depth_type: str = "D90",
                 exclude_aus: bool = True, exclude_zero: bool = True, sort_sr: bool = False,
                 segment_distance_tolerance: float = 100., smoothing_n_refinements: int = 5):

        self._smoothing_n = smoothing_n_refinements
        super(LeapfrogMultiFault, self).__init__(fault_geodataframe=fault_geodataframe, 
                                                 exclude_region_polygons=exclude_region_polygons,
                                                 exclude_region_min_sr=exclude_region_min_sr,
                                                 include_names=include_names,
                                                 depth_type=depth_type,
                                                 exclude_aus=exclude_aus, exclude_zero=exclude_zero, sort_sr=sort_sr)

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

    def add_fault(self, series: pd.Series, depth_type: str = "D90"):
        cfmFault = LeapfrogFault.from_series(series, parent_multifault=self, depth_type=depth_type)
        cfmFault.smoothed_trace = smooth_trace(cfmFault.nztm_trace, n_refinements=self.smoothing_n)
        self.faults.append(cfmFault)

    @property
    def smoothing_n(self):
        return self._smoothing_n

    @property
    def segment_distance_tolerance(self):
        return self._segment_distance_tolerance

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
                name_ls.append(fault.overall_name)
                sr_ls.append(max([seg.sr_best for seg in fault.segments]))

        df = pd.DataFrame({"name": name_ls, "sr": sr_ls})
        df.sort_values(["sr", "name"], ascending=(False, True), inplace=True)

        out_names = list(df.name)
        out_file_name = prefix + "_suggested_hierarchy.csv"

        with open(out_file_name, "w") as out_id:
            for name in out_names:
                out_id.write(name + "\n")





    @property
    def names(self):
        return [fault.name for fault in self.curated_faults]

    def find_connections(self):
        connections = []
        neighbour_connections = []
        for fault in self.faults:
            for other_fault in self.faults:
                if other_fault.name != fault.name:
                    if fault.nztm_trace.distance(other_fault.nztm_trace) <= self.segment_distance_tolerance:
                        print(f"Connection: {fault.name} and {other_fault.name}")
                        connections.append([fault.name, other_fault.name])
                        conditions = []
                        for p1, p2 in product([fault.end1, fault.end2], [other_fault.end1, other_fault.end2]):
                            conditions.append(p1.distance(p2) <= self.segment_distance_tolerance)
                        if any(conditions):
                            neighbour_connections.append([fault.name, other_fault.name])

        self._connections = connections
        self._neighbour_connections = neighbour_connections

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
                segs = [element.strip() for element in elements[1:]]
                cfault = ConnectedFaultSystem(overall_name=name, cfm_faults=self, segment_names=segs)
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
    def from_shp(cls, filename: str, exclude_region_polygons: List[Polygon] = None, depth_type: str = "D90",
                 exclude_region_min_sr: float = 1.8, sort_sr: bool = False, smoothing_n_refinements: int = 5):
        assert os.path.exists(filename)
        fault_geodataframe = gpd.GeoDataFrame.from_file(filename)
        multi_fault = cls(fault_geodataframe, exclude_region_polygons=exclude_region_polygons,
                          exclude_region_min_sr=exclude_region_min_sr, depth_type=depth_type, sort_sr=sort_sr,
                          smoothing_n_refinements=smoothing_n_refinements)
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


class LeapfrogFault(CfmFault):
    def __init__(self, parent_multifault: LeapfrogMultiFault = None, smoothing: int = 5,
                 trimming_gradient: float = 1.0, segment_distance_tolerance: float = 100.):
        self._end1 = None
        self._end2 = None
        self._neighbouring_segments = None
        self._neighbour_dict = {}
        self._neighbour_angle_dict = {}

        super(LeapfrogFault, self).__init__(None)
        self._parent = parent_multifault

        self.connections = []
        self._is_segment = False
        self._smoothing = smoothing
        self._trimming_gradient = trimming_gradient
        self._segment_distance_tolerance = segment_distance_tolerance
        self._contours = None
        self._smoothed_trace = None
        self._footprint = None

    @property
    def is_segment(self):
        return self._is_segment

    @property
    def smoothing(self):
        return self._smoothing

    @property
    def trimming_gradient(self):
        return self._trimming_gradient

    @property
    def neighbouring_segments(self):
        return self._neighbouring_segments

    @neighbouring_segments.setter
    def neighbouring_segments(self, segment_list: list):
        assert isinstance(segment_list, list)
        assert len(segment_list) <= 2
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
        return self._parent

    def depth_contour(self, depth: float, smoothing: bool = True, damping: int = None, km= False,
                      distance_tolerance: float = 100.):
        if depth <= 0:
            shift = depth / self.down_dip_vector[-1]
        else:
            shift = (-1 * depth) / self.down_dip_vector[-1]
        if km:
            shift *= 1000.

        xo, yo, zo = shift * self.down_dip_vector

        # if damping is not None:
        #     damped_contour = straighten(contour, self.dip_dir, damping)
        # else:
        #     damped_contour = contour

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
        contours = [self.depth_contour(depth, smoothing, damping, km) for depth in depths]

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
            trace = list(trace)[0]

        new_trace = LineString([(xi, yi, 0.) for xi, yi in trace.coords])

        self._nztm_trace = new_trace
        new_coord_array = np.array(new_trace)
        self._end1 = Point(new_coord_array[0])
        self._end2 = Point(new_coord_array[-1])


    @property
    def end1(self):
        return self._end1

    @property
    def end2(self):
        return self._end2

    def clipping_box(self, centre_point: Point, along_half_width: float, across_half_width: float = 100000.):
        point_array = np.array(centre_point)
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
        return self._segment_distance_tolerance

    @property
    def footprint(self):
        if self._footprint is None:
            self.calculate_footprint()
        return self._footprint

    def calculate_footprint(self, smoothed: bool = True, buffer: float = 15000.):
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

    def adjust_footprint(self, line_length: float = 1.5e5):
        terms = list(set(chain(*self.find_terminations())))
        if terms:
            cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
            fp_to_merge = [self.footprint] + [fault.footprint for fault in cutting_faults]
            merged_footprints = unary_union(fp_to_merge)
            cutting_lines = []
            for end_i, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
                if not any([end_i.distance(fault.nztm_trace) < self.segment_distance_tolerance for fault in cutting_faults]):
                    cutting_lines.append((LineString([np.array(end_i) + self.across_strike_vector * line_length,
                                                      np.array(end_i) - self.across_strike_vector * line_length]),
                                          other_end))
            if len(cutting_lines):
                if len(cutting_lines) > 1:
                    print(f"{self.name}: more than one cutting line, choosing first...")
                splitter, other_end = cutting_lines[0]
                split_footprint = split(merged_footprints, splitter)
                kept_polys = [poly for poly in list(split_footprint) if other_end.within(poly)]
                if len(kept_polys) > 1:
                    print(f"{self.name}: more than one cut polygon, choosing first...")
                self._footprint = kept_polys[0]

            else:
                self._footprint = merged_footprints

    def extend_footprint(self, end_i: Point, other_end: Point, distance: float = 40.e3):
        diff_vector = np.array(other_end) - np.array(end_i)
        if np.dot(diff_vector, self.along_strike_vector) > 0:
            strike_direction = -1 * self.along_strike_vector
        else:
            strike_direction = self.along_strike_vector

        bottom_trace = list(self.contours.geometry)[-1]














class LeapfrogConnectedFault:
    def __init__(self, fault_list: List[LeapfrogFault], parent_multifault: LeapfrogMultiFault = None):
        assert all([isinstance(x, LeapfrogFault) for x in fault_list])
        self._parent = parent_multifault
        self.faults = fault_list

        self.connections = []




    def join_combine_surface_traces(self, smoothing: int = 0, simplification: int = 0):
        pass

    def depth_contours(self, smoothing: int = 0, simplification: int = 0, angle_gap: Union[float, int] = None):
        pass


class LeapfrogFaultSegment(LeapfrogFault):
    def __init__(self, parent_fault: LeapfrogConnectedFault, parent_multifault: LeapfrogMultiFault = None, ):
        super(LeapfrogFault, self).__init__(parent_multifault)
        self._fault = parent_fault
        self.neighbouring_segments = []




class LeapfrogFaultModel:
    def __init__(self, fault_list: List[Union[LeapfrogFault, LeapfrogConnectedFault]]):
        assert all([isinstance(x, (LeapfrogFault, LeapfrogConnectedFault)) for x in fault_list])
        self.faults = fault_list
