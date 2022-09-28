"""
Module for connecting segments into a connected fault system
"""
from typing import Union, List
import fnmatch
from itertools import chain

import numpy as np
from shapely.geometry import MultiLineString, LineString, Point, Polygon
from shapely.ops import linemerge, unary_union
import geopandas as gpd

from eq_fault_geom.geomio.cfm_faults import smallest_difference
from fault_mesh.utilities.smoothing import smooth_trace, merge_multiple_nearly_adjacent_segments, align_two_nearly_adjacent_segments
from fault_mesh.utilities.cutting import cut_line_at_multiple_points, cut_line_at_point


class ConnectedFaultSystem:
    def __init__(self, overall_name: str, cfm_faults, segment_names: list = None,
                 search_patterns: Union[str, list] = None,
                 excluded_names: Union[str, list] = None, tolerance: float = 100.,
                 smooth_trace_refinements: int = 5):
        self.name = overall_name
        self._overall_trace = None
        self._contours = None
        self._footprint = None
        segments = []

        assert any([segment_names is not None, search_patterns is not None])

        if segment_names is not None:
            assert isinstance(segment_names, list)
            assert len(segment_names)
            segment_names_list = segment_names
        else:
            segment_names_list = []

        if search_patterns is not None:
            if isinstance(search_patterns, str):
                search_pattern_list = [search_patterns]
            else:
                assert isinstance(search_patterns, list)
                assert len(search_patterns)
                search_pattern_list = search_patterns

        if isinstance(excluded_names, str):
            excluded_list = [excluded_names]
        elif excluded_names is None:
            excluded_list = []
        else:
            assert isinstance(excluded_names, list)
            assert len(excluded_names)
            excluded_list = excluded_names

        if len(excluded_list):
            assert not any([x in segment_names_list for x in excluded_list])

        for fault in cfm_faults.faults:
            name = fault.name
            if search_patterns is not None:
                if any([fnmatch.fnmatch(name, pattern) for pattern in search_pattern_list]):
                    if not any([fnmatch.fnmatch(name, pattern) for pattern in excluded_list]):
                        segments.append(fault)

            if segment_names is not None:
                if any([fnmatch.fnmatch(name, pattern) for pattern in segment_names_list]):
                    segments.append(fault)

        for segment in segments:
            other_segments = MultiLineString([other_seg.nztm_trace for other_seg in segments if other_seg != segment])
            closest_dist = segment.nztm_trace.distance(other_segments)
            if closest_dist > tolerance:
                raise ValueError(f"Fault traces >{tolerance} m apart: {segment.name}")
            else:
                neighbour_list = []
                for other_seg in segments:
                    if other_seg != segment:
                        if segment.nztm_trace.distance(other_seg.nztm_trace) < tolerance:
                            neighbour_list.append(other_seg)

                segment.neighbouring_segments = neighbour_list
            segment.parent_connected_fault = self

        # Order segments
        # Find furthest south or west

        sorted_segments_s_to_n = sorted(segments, key=lambda x: (x.nztm_trace.bounds[0],
                                                                 x.nztm_trace.bounds[1]))
        sorted_segments_w_to_e = sorted(segments, key=lambda x: (x.nztm_trace.bounds[1],
                                                                 x.nztm_trace.bounds[0]))

        if sorted_segments_s_to_n[0] != sorted_segments_w_to_e[0]:
            print(f"Warning: southmost and westernmost segments differ: {overall_name}")
            print(f"Setting {sorted_segments_s_to_n[0].name} as first segment")
        first_segment = sorted_segments_s_to_n[0]

        ordered_segments = [first_segment]
        while len(ordered_segments) < len(segments):
            for segment in segments:
                if segment not in ordered_segments:
                    if segment.nztm_trace.distance(ordered_segments[-1].nztm_trace) <= tolerance:
                        ordered_segments.append(segment)
        self.segments = ordered_segments
        self.overall_trace = linemerge([seg.nztm_trace for seg in self.segments])

        # boundaries = [l1.nztm_trace.intersection(l2.nztm_trace) for l1, l2 in zip(self.segments, self.segments[1:])]
        boundaries = []
        for l1, l2 in zip(self.segments, self.segments[1:]):
            if l1.nztm_trace.distance(l2.nztm_trace) == 0.:
                boundaries.append(l1.nztm_trace.intersection(l2.nztm_trace))
            else:
                new_l1, new_l2 = align_two_nearly_adjacent_segments([l1.nztm_trace, l2.nztm_trace])
                boundaries.append(new_l1.intersection(new_l2))



        self.segment_boundaries = [bound for bound in boundaries if not bound.is_empty]


        if smooth_trace_refinements is not None:
            self.smoothed_overall_trace = smooth_trace(self.overall_trace, n_refinements=smooth_trace_refinements)
            if len(self.segment_boundaries) >= 2:
                self.smoothed_segments = cut_line_at_multiple_points(self.smoothed_overall_trace, self.segment_boundaries)
            else:
                self.smoothed_segments = cut_line_at_point(self.smoothed_overall_trace, self.segment_boundaries[0])
            for smoothed_segment in self.smoothed_segments:
                closest_segment = min(self.segments, key=lambda x: x.nztm_trace.centroid.distance(smoothed_segment.centroid))
                closest_segment.smoothed_trace = smoothed_segment
        else:
            self.smoothed_overall_trace = None
            self.smoothed_segments = None

        self.parent = self.segments[0].parent

    @property
    def segment_names(self):
        return [segment.name for segment in self.segments]

    @property
    def trace(self):
        return self.overall_trace

    @property
    def nztm_trace(self):
        return self.overall_trace

    @property
    def smoothed_trace(self):
        return self.smoothed_overall_trace

    def depth_contour(self, depth: float, smoothing: bool = True, damping: int = None, km: bool = False):
        contours = [segment.depth_contour(depth, smoothing) for segment in self.segments]
        valid_contours = [contour for contour in contours if contour is not None]
        return MultiLineString(valid_contours)

    def generate_depth_contours(self, depths: Union[np.ndarray, List[float]], smoothing: bool = True, damping: int = None,
                       km: bool = False):
        contours = [self.depth_contour(depth, smoothing) for depth in depths]

        if max(depths) > 0:
            depths *= -1

        self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours)

    @property
    def contours(self):
        return self._contours

    @contours.setter
    def contours(self, contours: gpd.GeoDataFrame):
        self._contours = contours



        # Might need to add interpolation between nearby but not identical segment ends. Do this when we get to it.

    @property
    def overall_trace(self):
        return self._overall_trace

    @overall_trace.setter
    def overall_trace(self, trace: Union[LineString, MultiLineString]):
        if isinstance(trace, LineString):
            self._overall_trace = trace
        else:
            assert isinstance(trace, MultiLineString)
            self._overall_trace = merge_multiple_nearly_adjacent_segments(list(trace.geoms))

    def trace_and_contours(self, smoothed: bool = True):
        assert self.contours is not None
        if smoothed:
            return unary_union([self.smoothed_trace] + list(self.contours.geometry) + list(self.end_lines(smoothed=smoothed).geoms))
        else:
            return unary_union([self.overall_trace] + list(self.contours.geometry) + list(self.end_lines(smoothed=smoothed).geoms))

    def end_polygon(self, smoothed: bool = False, distance: float= 1.e5):
        if smoothed:
            end1 = Point(self.smoothed_overall_trace.coords[0])
            end2 = Point(self.smoothed_overall_trace.coords[-1])
        else:
            end1 = Point(self.trace.coords[0])
            end2 = Point(self.trace.coords[-1])

        seg1 = min(self.segments, key=lambda x: x.nztm_trace.centroid.distance(end1))
        seg2 = min(self.segments, key=lambda x: x.nztm_trace.centroid.distance(end2))

        if smoothed:
            assert end1.distance(seg1.smoothed_trace) < 500.
            assert end2.distance(seg2.smoothed_trace) < 500.
        else:
            assert end1.distance(seg1.nztm_trace) < 500.
            assert end2.distance(seg2.nztm_trace) < 500.

        line1 = LineString(np.vstack([np.array(end1.coords) + distance * seg1.across_strike_vector,
                                      np.array(end1.coords) - distance * seg1.across_strike_vector]))

        if smallest_difference(seg1.dip_dir, seg2.dip_dir) > 90:
            line2 = LineString(np.vstack([np.array(end2.coords) + distance * seg2.across_strike_vector,
                                          np.array(end2.coords) - distance * seg2.across_strike_vector]))
        else:
            line2 = LineString(np.vstack([np.array(end2.coords) - distance * seg2.across_strike_vector,
                                          np.array(end2.coords) + distance * seg2.across_strike_vector]))

        if line1.intersects(line2):
            intersection = line1.intersection(line2)
            p1 = max([p for p in line1.coords], key=lambda x: Point(x).distance(intersection))
            p2 = max([p for p in line2.coords], key=lambda x: Point(x).distance(intersection))
            edge_poly = Polygon([p1, intersection, p2])

        else:
            edge_poly = Polygon(np.vstack([np.array(line1.coords), np.array(line2.coords)]))

        return edge_poly

    @property
    def footprint(self):
        if self._footprint is None:
            self.calculate_footprint()
        return self._footprint

    def calculate_footprint(self, smoothed: bool = True, buffer: float = 15000.):
        footprint = self.trace_and_contours(smoothed).minimum_rotated_rectangle.buffer(buffer, cap_style=2).intersection(self.end_polygon(smoothed=smoothed))
        self._footprint = footprint

    def end_lines(self, smoothed: bool = False, depth: float = 20.e3, spacing: float = 2000.):
        if smoothed:
            end1 = Point(np.array(self.smoothed_overall_trace.coords[0]))
            end2 = Point(np.array(self.smoothed_overall_trace.coords[-1]))
        else:
            end1 = Point(np.array(self.trace.coords[0]))
            end2 = Point(np.array(self.trace.coords[-1]))

        seg1 = min(self.segments, key=lambda x: x.nztm_trace.centroid.distance(end1))
        seg2 = min(self.segments, key=lambda x: x.nztm_trace.centroid.distance(end2))

        end_line_list = []
        for end, seg in zip([end1, end2], [seg1, seg2]):
            vertex_list = [np.array(end.coords)]
            distance = depth / np.sin(np.radians(seg.dip_best))
            for dist in np.arange(spacing, distance, spacing):
                vertex_list.append(np.array(end.coords) + dist * seg.down_dip_vector)

            vertex_list.append(np.array(end.coords) + distance * seg.down_dip_vector)
            end_line_list.append(LineString(np.vstack(vertex_list)))

        combined = unary_union(end_line_list)
        if isinstance(combined, LineString):
            return MultiLineString([combined])
        else:
            return combined

    def find_terminations(self):
        return self.parent.find_terminations(self.name)


    def adjust_footprint(self, smoothed: bool = True):
        if smoothed:
            end1 = Point(self.smoothed_overall_trace.coords[0])
            end2 = Point(self.smoothed_overall_trace.coords[-1])
        else:
            end1 = Point(self.trace.coords[0])
            end2 = Point(self.trace.coords[-1])
        terms = list(set(chain(*self.find_terminations())))
        if terms:
            cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
            for fault in cutting_faults:
                for nearest_end, other_end in zip([end1, end2], [end2, end1]):
                    nearest_seg_this_fault = min(self.segments, key=lambda x: x.nztm_trace.distance(nearest_end))
                    if nearest_end.distance(fault.nztm_trace) < 1.e3:
                        if isinstance(fault, ConnectedFaultSystem):
                            closest_seg = min(fault.segments, key=lambda x: x.nztm_trace.distance(nearest_end))
                        else:
                            closest_seg = fault

                        nearest_seg_this_fault.extend_footprint(nearest_end, other_end, closest_seg)
    # def adjust_footprint(self, line_length: float = 1.5e5, smoothed: bool = True):
    #     terms = list(set(chain(*self.find_terminations())))
    #     if terms:
    #         cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
    #         fp_to_merge = [self.footprint] + [fault.footprint for fault in cutting_faults]
    #         merged_footprints = unary_union(fp_to_merge)
    #         cutting_lines = []
    #
    #         if smoothed:
    #             end1 = Point(self.smoothed_overall_trace.coords[0])
    #             end2 = Point(self.smoothed_overall_trace.coords[-1])
    #         else:
    #             end1 = Point(self.trace.coords[0])
    #             end2 = Point(self.trace.coords[-1])
    #
    #         for end_i, other_end in zip([end1, end2], [end2, end1]):
    #             nearest_seg = min(self.segments, key=lambda x:x.nztm_trace.centroid.distance(end_i))
    #             if not any([end_i.distance(fault.nztm_trace) < self.parent.segment_distance_tolerance for fault in cutting_faults]):
    #                 cutting_lines.append((LineString([np.array(end_i) + nearest_seg.across_strike_vector * line_length,
    #                                                   np.array(end_i) - nearest_seg.across_strike_vector * line_length]),
    #                                       other_end))
    #         if len(cutting_lines):
    #             if len(cutting_lines) > 1:
    #                 print(f"{self.name}: more than one cutting line, choosing first...")
    #             splitter, other_end = cutting_lines[0]
    #             split_footprint = split(merged_footprints, splitter)
    #             kept_polys = [poly for poly in list(split_footprint) if other_end.within(poly)]
    #             if len(kept_polys) > 1:
    #                 print(f"{self.name}: more than one cut polygon, choosing first...")
    #             self._footprint = kept_polys[0]
    #
    #         else:
    #             self._footprint = merged_footprints





