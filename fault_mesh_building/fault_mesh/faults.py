from __future__ import annotations
import os
from typing import List, Union
from itertools import product

import numpy as np
from eq_fault_geom.geomio.cfm_faults import CfmMultiFault, CfmFault, normalize_bearing, smallest_difference
import geopandas as gpd
import pandas as pd

from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString, Point, Polygon

from fault_mesh.smoothing import straighten, smooth_trace
from fault_mesh.utilities.cutting import cut_line_between_two_points, cut_line_at_point
from fault_mesh.utilities.graph import connected_nodes, suggest_combined_name


class LeapfrogMultiFault(CfmMultiFault):
    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, exclude_region_polygons: list = None,
                 exclude_region_min_sr: float = 1.8, include_names: list = None, depth_type: str = "D90",
                 exclude_aus: bool = True, exclude_zero: bool = True, sort_sr: bool = False,
                 segment_distance_tolerance: float = 100.):
        super(LeapfrogMultiFault, self).__init__(fault_geodataframe=fault_geodataframe, 
                                                 exclude_region_polygons=exclude_region_polygons,
                                                 exclude_region_min_sr=exclude_region_min_sr,
                                                 include_names=include_names,
                                                 depth_type=depth_type,
                                                 exclude_aus=exclude_aus, exclude_zero=exclude_zero, sort_sr=sort_sr)

        self._segment_distance_tolerance = segment_distance_tolerance
        self._cutting_hierarchy = []
        self._includes_connected = False
        self._curated_faults = None
        self._connections = None
        self._neighbour_connections = None

    def add_fault(self, series: pd.Series, depth_type: str = "D90"):
        cfmFault = LeapfrogFault.from_series(series, parent_multifault=self, depth_type=depth_type)
        self.faults.append(cfmFault)

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
    def curated_faults(self):
        if self._curated_faults is None:
            return self.faults
        else:
            return self._curated_faults

    def read_cutting_hierarchy(self, hierarchy_file: str):
        assert os.path.exists(hierarchy_file)
        with open(hierarchy_file, "r") as hfile:
            names = hfile.readlines()
            self.cutting_hierarchy = list(names)

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

    def suggest_fault_systems(self, out_prefix):
        connected = connected_nodes(self.neighbour_connections)
        for connected_set in connected:
            suggested_name = suggest_combined_name(connected_set)
            out_list = [suggested_name] + list(connected_set)
            out_str = ",".join(out_list) + "\n"





class LeapfrogFault(CfmFault):
    def __init__(self, parent_multifault: LeapfrogMultiFault = None, smoothing: int = None,
                 trimming_gradient: float = None, segment_distance_tolerance: float = 100.):
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
                self._neighbour_dict[self.end1] = seg
                self._neighbour_angle_dict[self.end1] = self.neighbour_angle(seg)
            elif seg.nztm_trace.distance(self.end2) <= self.segment_distance_tolerance:
                self._neighbour_dict[self.end2] = seg
                self._neighbour_angle_dict[self.end2] = self.neighbour_angle(seg)
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

    def depth_contour(self, depth: float, smoothing: int = None, damping: int = None, km= False):
        if depth <= 0:
            shift = depth / self.down_dip_vector[-1]
        else:
            shift = (-1 * depth) / self.down_dip_vector[-1]
        if km:
            shift *= 1000.

        xo, yo, zo = shift * self.down_dip_vector
        contour = translate(self.nztm_trace, xo, yo, zo)

        if damping is not None:
            damped_contour = straighten(contour, self.dip_dir, damping)
        else:
            damped_contour = contour

        if smoothing is not None:
            smoothed_contour = smooth_trace(damped_contour, smoothing)

        elif self.smoothing is not None:
            smoothed_contour = smooth_trace(damped_contour, self.smoothing)

        else:
            smoothed_contour = damped_contour

        if self.is_segment:
            if len(self.neighbouring_segments) > 1:
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                e1_box_intersection = e1_box.intersection(self.nztm_trace)
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                e2_box_intersection = e2_box.intersection(self.nztm_trace)

                trimmed_contour = cut_line_between_two_points(smoothed_contour, [e1_box_intersection,
                                                                                 e2_box_intersection])

            elif self.end1 in self.neighbour_dict.keys():
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                e1_box_intersection = e1_box.intersection(self.nztm_trace)
                split_line = cut_line_at_point(smoothed_contour, e1_box_intersection)

                trimmed_contour = split_line[1]

            else:  # self.end2 in self.neighbour_dict.keys()
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                e2_box_intersection = e2_box.intersection(self.nztm_trace)
                split_line = cut_line_at_point(smoothed_contour, e2_box_intersection)

                trimmed_contour = split_line[0]

        else:
            trimmed_contour = smoothed_contour

        return trimmed_contour

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

    def clipping_box(self, centre_point: Point, along_half_width: float, across_half_width: float = 10000.):
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
        end_angle = self.neighbour_angle_dict[end_i]
        end_width = gradient_adjustment * np.tan(np.radians(end_angle)) / (depth / np.sin(np.radians(end_angle)))

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
        return max([dip_diff, strike_diff]) / 2.





    @property
    def segment_distance_tolerance(self):
        return self._segment_distance_tolerance









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
