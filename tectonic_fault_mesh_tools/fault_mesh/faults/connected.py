"""
Module for connecting segments into a connected fault system
"""
from typing import Union, List
import fnmatch
from itertools import chain

from build.lib.fault_mesh.utilities import smoothing
import numpy as np
from shapely.geometry import MultiLineString, LineString, Point, Polygon
from shapely.ops import linemerge, unary_union
from shapely import get_point
from shapely.affinity import translate
import geopandas as gpd
import matplotlib.pyplot as plt
import triangle as tr
import meshio
from scipy.interpolate import RBFInterpolator

from fault_mesh.faults.generic import smallest_difference, calculate_dip_direction, reverse_bearing, dip_direction_ranges, reverse_line
from fault_mesh.utilities.smoothing import smooth_trace
from fault_mesh.utilities.merging import merge_multiple_nearly_adjacent_segments, align_two_nearly_adjacent_segments
from fault_mesh.utilities.cutting import cut_line_at_multiple_points, cut_line_at_point
from fault_mesh.utilities.meshing import get_strike_dip_from_normal, fit_plane_to_points, weighted_circular_mean, most_common_or_first, triangulate_contours
from fault_mesh.utilities.splines import spline_fit_contours


class ConnectedFaultSystem:
    """_summary_
    """
    def __init__(self, overall_name: str, cfm_faults, segment_names: list = None,
                 search_patterns: Union[str, list] = None,
                 excluded_names: Union[str, list] = None, tolerance: float = 100.,
                 smooth_trace_refinements: int = 5, trimming_gradient: float = 1.):
        """A class representing a connected fault system composed of multiple fault segments.
        
        This class handles the integration of individual fault segments into a connected system,
        including segment alignment, boundary identification, and trace smoothing.

        :param overall_name: Name identifier for the connected fault system
        :type overall_name: str
        :param cfm_faults: Collection of fault objects to extract segments from
        :type cfm_faults: FaultCollection
        :param segment_names: Explicit list of segment names to include, defaults to None
        :type segment_names: list, optional
        :param search_patterns: Pattern(s) to match segment names against, defaults to None
        :type search_patterns: Union[str, list], optional
        :param excluded_names: Segment names to explicitly exclude, defaults to None
        :type excluded_names: Union[str, list], optional
        :param tolerance: Maximum distance between segments to be considered connected in meters, defaults to 100.
        :type tolerance: float, optional
        :param smooth_trace_refinements: Number of iterations for trace smoothing algorithm, defaults to 5
        :type smooth_trace_refinements: int, optional
        :param trimming_gradient: Gradient used for trimming fault segments, defaults to 1.
        :type trimming_gradient: float, optional
        :raises ValueError: If fault segments are further apart than the specified tolerance
        """
        self.name = overall_name
        self._overall_trace = None
        self._contours = None
        self._footprint = None
        self._segment_distance_tolerance = tolerance
        self._trimming_gradient = trimming_gradient
        self._smooth_trace_refinements = smooth_trace_refinements
        self._dip_dir = None
        segments = []

        assert any([segment_names is not None, search_patterns is not None])

        # Handle segment_names parameter
        if segment_names is not None:
            assert isinstance(segment_names, list)
            assert len(segment_names)  # Ensure the list is not empty
            segment_names_list = segment_names
        else:
            segment_names_list = []

        # Handle search_patterns parameter
        if search_patterns is not None:
            if isinstance(search_patterns, str):
                search_pattern_list = [search_patterns]  # Convert single string pattern to list
            else:
                assert isinstance(search_patterns, list)
                assert len(search_patterns)  # Ensure the list is not empty
                search_pattern_list = search_patterns

        # Handle excluded_names parameter
        if isinstance(excluded_names, str):
            excluded_list = [excluded_names]  # Convert single string exclusion to list
        elif excluded_names is None:
            excluded_list = []  # Initialize empty list if no exclusions
        else:
            assert isinstance(excluded_names, list)
            assert len(excluded_names)  # Ensure the list is not empty
            excluded_list = excluded_names

        # Verify no overlap between segment_names and excluded_names
        if len(excluded_list):
            assert not any([x in segment_names_list for x in excluded_list])

        # Collect fault segments based on search patterns and segment names
        for fault in cfm_faults.faults:
            name = fault.name
            if search_patterns is not None:
            # Check if fault name matches any of the search patterns
                if any([fnmatch.fnmatch(name, pattern) for pattern in search_pattern_list]):
                    # Ensure fault name is not in the exclusion list
                    if not any([fnmatch.fnmatch(name, pattern) for pattern in excluded_list]):
                        segments.append(fault)

            if segment_names is not None:
            # Add fault if its name matches any in the segment_names list
                if any([fnmatch.fnmatch(name, pattern) for pattern in segment_names_list]):
                    segments.append(fault)

        # Verify segment connectivity and store neighbor information
        for segment in segments:
            # Create a MultiLineString of all other segments' traces
            other_segments = MultiLineString([other_seg.nztm_trace for other_seg in segments if other_seg != segment])
            closest_dist = segment.nztm_trace.distance(other_segments)
            
            # Check if segment is within tolerance distance of at least one other segment
            if closest_dist > tolerance:
                raise ValueError(f"Fault traces >{tolerance} m apart: {segment.name}")
            else:
            # Identify neighboring segments within tolerance distance
                neighbour_list = []
                for other_seg in segments:
                    if other_seg != segment:
                        if segment.nztm_trace.distance(other_seg.nztm_trace) < tolerance:
                            neighbour_list.append(other_seg)
                segment.neighbouring_segments = neighbour_list
                
            # Set reference back to this connected fault system
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
        if self._smooth_trace_refinements is None:
            self.overall_trace = merge_multiple_nearly_adjacent_segments([seg.nztm_trace for seg in self.segments],
                                                                         densify=None)
        else:
            self.overall_trace = merge_multiple_nearly_adjacent_segments([seg.nztm_trace for seg in self.segments])

        # boundaries = [l1.nztm_trace.intersection(l2.nztm_trace) for l1, l2 in zip(self.segments, self.segments[1:])]
        boundaries = []
        for l1, l2 in zip(self.segments, self.segments[1:]):
            if l1.nztm_trace.distance(l2.nztm_trace) == 0.:
                boundaries.append(l1.nztm_trace.intersection(l2.nztm_trace))
            else:
                if self._smooth_trace_refinements is None:
                    new_l1, new_l2 = align_two_nearly_adjacent_segments([l1.nztm_trace, l2.nztm_trace], tolerance=tolerance,
                                                                        densify=None)
                else:
                    new_l1, new_l2 = align_two_nearly_adjacent_segments([l1.nztm_trace, l2.nztm_trace],
                                                                        tolerance=tolerance)
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
    def segment_distance_tolerance(self):
        return self._segment_distance_tolerance

    @property
    def trace(self):
        return self.overall_trace

    @property
    def nztm_trace(self):
        return self.overall_trace

    @property
    def nztm_trace_geoseries(self):
        if self.parent.epsg is not None:
            return gpd.GeoSeries(self.nztm_trace, crs=self.parent.epsg)
        else:
            return gpd.GeoSeries(self.nztm_trace)

    @property
    def smoothed_trace(self):
        return self.smoothed_overall_trace
    
    @property
    def dip_dir(self):
        return self._dip_dir

    def depth_contour(self, depth: float, smoothing: bool = True, damping: int = None, km: bool = False):
        contours = [segment.depth_contour(depth, smoothing) for segment in self.segments]
        valid_contours = [contour for contour in contours if contour is not None]
        return MultiLineString(valid_contours)

    def generate_depth_contours(self, depths: Union[np.ndarray, List[float]], smoothing: bool = True, damping: int = None,
                       km: bool = False):
        contours = [self.depth_contour(depth, smoothing) for depth in depths]

        if max(depths) > 0:
            depths *= -1
        if self.parent.epsg is not None:
            self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours, crs=self.parent.epsg)
        else:
            self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours)

    def generate_simple_contours(self, depths: Union[np.ndarray, List[float]], km: bool = False):
        contours = [self.depth_contour(depth) for depth in depths]
        seg_lengths = [seg.nztm_trace.length for seg in self.segments]
        seg_dips = [seg.dip_best for seg in self.segments]
        average_dip = weighted_circular_mean(seg_dips, seg_lengths)
        self.calculate_overall_dip_direction()

        if self.dip_dir is None:
            # Assume vertical
            dd_vec = np.array([0., 0., -1])
        else:
            z = np.sin(np.radians(average_dip))
            x, y = np.cos(np.radians(average_dip)) * np.array([np.sin(np.radians(self.dip_dir)),
                                                                np.cos(np.radians(self.dip_dir))])
            dd_vec = np.array([x, y, -z])

        contours = []
        for depth in depths:
            if depth <= 0:
                shift = depth / dd_vec[-1]
            else:
                shift = (-1 * depth) / dd_vec[-1]
            if km:
                shift *= 1000.

            xo, yo, zo = shift * dd_vec

            shifted_contour = translate(self.nztm_trace, xo, yo, zo)
            contours.append(shifted_contour)

        if max(depths) > 0:
            depths *= -1
        if self.parent.epsg is not None:
            simple_contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours, crs=self.parent.epsg)
        else:
            simple_contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours)

        return simple_contours

    def calculate_overall_dip_direction(self, tolerance: float = 10.):
        """
        Compares dip direction string (e.g. NW) with
        :return:
        """
        dip_dir_strings = [seg.dip_dir_str for seg in self.segments if seg.dip_dir_str is not None]
        dd_str = most_common_or_first(dip_dir_strings)
        if any([a is None for a in [dd_str, self.nztm_trace]]):
            print("Insufficient information to validate dip direction")
            if self.nztm_trace is not None and dd_str is None:
                dd_from_trace = calculate_dip_direction(self.nztm_trace)
            self._dip_dir = dd_from_trace
        else:
            # Trace and dip direction

            dd_from_trace = calculate_dip_direction(self.nztm_trace)
            if dd_str != "SUBVERTICAL AND VARIABLE":
                min_dd_range, max_dd_range = dip_direction_ranges[dd_str]
                if dd_str != "N":
                    if not all([min_dd_range - tolerance <= dd_from_trace, dd_from_trace <= max_dd_range + tolerance]):
                        reversed_dd = reverse_bearing(dd_from_trace)
                        if all([min_dd_range - tolerance <= reversed_dd, reversed_dd <= max_dd_range + tolerance]):
                            self._nztm_trace = reverse_line(self.nztm_trace)
                            self._dip_dir = reversed_dd
                        else:
                            print("{}: Supplied trace and dip direction {} are inconsistent: expect either {:.1f}"
                                  "or {:.1f} dip azimuth. Please check...".format(self.name, self.dip_dir_str,
                                                                                  dd_from_trace, reversed_dd))

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
    
    def generate_sr_rake_points(self, depths: Union[np.ndarray, List[float]], smoothing: bool = True,
                                km: bool = False):
        
        """
        Generate strike and rake points for the fault system at specified depths.
        
        :param depths: Depths at which to generate strike and rake points
        :type depths: Union[np.ndarray, List[float]]
        :param smoothing: Whether to apply smoothing to the traces, defaults to True
        :type smoothing: bool, optional
        :param km: Whether depths are in kilometers, defaults to False
        :type km: bool, optional
        """
        point_geoms = []
        point_depths = []
        point_srs = []
        point_rakes = []
        for depth in depths:
            if km:
                depth *= 1000.
            for segment in self.segments:
                contour = segment.depth_contour(depth, smoothing=smoothing, km=km)
                if contour is not None:
                    for point_i in [0, -1]:
                        end_point = get_point(contour, point_i)
                        point_geoms.append(end_point)
                        point_depths.append(depth)
                        point_srs.append(segment.sr_best)
                        point_rakes.append(segment.rake_best)
        if self.parent.epsg is not None:
            sr_rake_points = gpd.GeoDataFrame({
                "depth": point_depths,
                "slip_rate": point_srs,
                "rake": point_rakes
            }, geometry=point_geoms, crs=self.parent.epsg)  
        else:
            sr_rake_points = gpd.GeoDataFrame({
                "depth": point_depths,
                "slip_rate": point_srs,
                "rake": point_rakes
            }, geometry=point_geoms)

        return sr_rake_points

    @property
    def contours(self):
        return self._contours

    @contours.setter
    def contours(self, contours: gpd.GeoDataFrame):
        self._contours = contours

    @property
    def sampled_dip(self):
        """
        Dip value sampled from the fault trace. Used for calculating dip in depth contours.
        :return:
        """
        return self._sampled_dip
    
    @sampled_dip.setter
    def sampled_dip(self, dip: float):
        """
        Set the sampled dip value.
        :param dip: Dip value in degrees.
        """
        assert isinstance(dip, (int, float))
        assert 0 <= dip <= 1, "Dip must be between 0 and 1"
        self._sampled_dip = dip



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
            self._overall_trace = merge_multiple_nearly_adjacent_segments(list(trace.geoms),
                                                                          tolerance=self.segment_distance_tolerance)

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

    @property
    def footprint_geoseries(self):
        if self.parent.epsg is not None:
            return gpd.GeoSeries(self.footprint, crs=self.parent.epsg)
        else:
            return gpd.GeoSeries(self.footprint)

    @property
    def footprint_linestring(self):
        return LineString(self.footprint.exterior.coords)

    def calculate_footprint(self, buffer: float = 15000.):
        if self.parent.smoothing_n is not None:
            footprint = self.trace_and_contours(smoothed=True).minimum_rotated_rectangle.buffer(buffer, cap_style=2).intersection(self.end_polygon(smoothed=True))
        else:
            footprint = self.trace_and_contours(smoothed=False).minimum_rotated_rectangle.buffer(buffer, cap_style=2).intersection(self.end_polygon(smoothed=False))
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
            if not combined.is_empty:
                return MultiLineString([combined])
            else:
                assert not any(end_line.is_empty for end_line in end_line_list)
                return MultiLineString(end_line_list)
        else:
            return combined

    def find_terminations(self):
        return self.parent.find_terminations(self.name)

    def adjust_footprint(self):
        if self.parent.smoothing_n is not None:
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

    def adjust_trace(self, extend_distance: float = 20.e3, fit_distance: float = 5.e3, resolution: float = 1.e3):
        if self.parent.smoothing_n is not None:
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
                    if nearest_end.distance(fault.nztm_trace) < 1.e3:
                        nearest_seg_this_fault = min(self.segments, key=lambda x: x.nztm_trace.distance(nearest_end))
                        nearest_seg_this_fault.extend_trace(nearest_end, other_end, fit_distance=fit_distance, extend_distance=extend_distance, resolution=resolution)

        if self._smooth_trace_refinements is None:
            self.overall_trace = merge_multiple_nearly_adjacent_segments([seg.nztm_trace for seg in self.segments],
                                                                         densify=None)
        else:
            self.overall_trace = merge_multiple_nearly_adjacent_segments([seg.nztm_trace for seg in self.segments])


    def mesh_fault_surface(self, resolution: float = 1000., spline_resolution: float = 100., plane_fitting_eps: float = 1.0e-5, check_mesh: bool = False, check_strike_dip: bool = False):
        """Generate a triangular mesh representing the fault surface.
        
        This method creates a 3D triangular mesh of the fault surface by:
        1. Fitting splines to depth contours
        2. Fitting a plane to all contour points
        3. Projecting contours onto the plane
        4. Interpolating between contours using thin plate spline
        5. Triangulating the resulting surface
        
        :param resolution: Target spatial resolution of the mesh in meters, defaults to 1000.
        :type resolution: float, optional
        :param spline_resolution: Point spacing for spline fitting in meters, defaults to 100.
        :type spline_resolution: float, optional
        :param plane_fitting_eps: Tolerance for plane fitting algorithm, defaults to 1.0e-5
        :type plane_fitting_eps: float, optional
        :param check_mesh: Whether to display visualization of the mesh for debugging, defaults to False
        :type check_mesh: bool, optional
        :param check_strike_dip: Whether to print strike and dip information for debugging, defaults to False
        :type check_strike_dip: bool, optional
        :return: Triangular mesh representing the fault surface
        :rtype: meshio.Mesh
        """
        spline_contours = self.generate_simple_contours(self.contours.depth.values, km=False)
        spline_contour_list = [np.array(contour.coords) for contour in spline_contours.geometry.values]

        spline_contour_points = np.vstack(spline_contour_list)
        spline_contour_dict = {tuple(point): i for (i,point) in enumerate(spline_contour_points)}

        spline_contour_boundary = np.vstack([spline_contour_list[0], 
                                    np.vstack([spline_contour_list[i][-1] for i in range(1, len(spline_contour_list) -1)]),
                                    spline_contour_list[-1][::-1],
                                    np.vstack([spline_contour_list[i][0] for i in range(len(spline_contour_list) -1, 0, -1)])])

        contours_list = []
        for contour in self.contours.geometry.values:
            if isinstance(contour, LineString):
                contours_list.append(np.array(contour.coords))
            else:
                assert isinstance(contour, MultiLineString)
                contours_list.append(np.vstack([np.array(line.coords) for line in contour.geoms]))


        all_contour_points = np.vstack(contours_list)
        plane_normal, plane_origin = fit_plane_to_points(all_contour_points, eps=plane_fitting_eps)

        # Calculate strike and dip from normal
        strike, dip = get_strike_dip_from_normal(plane_normal)
        if check_strike_dip:
            print(f"Plane strike: {strike:.1f}°, dip: {dip:.1f}°")

        # Calculate in-plane x vector (along strike)
        strike_rad = np.radians(strike)
        in_plane_x_vector = np.array([np.sin(strike_rad), np.cos(strike_rad), 0])
        in_plane_x_vector = in_plane_x_vector - np.dot(in_plane_x_vector, plane_normal) * plane_normal
        in_plane_x_vector = in_plane_x_vector / np.linalg.norm(in_plane_x_vector)

        # # Calculate in-plane y vector (down dip)
        # in_plane_y_vector = np.cross(plane_normal, in_plane_x_vector)
        # in_plane_y_vector = in_plane_y_vector / np.linalg.norm(in_plane_y_vector)

        # set in-plane y vector to be straight down
        in_plane_y_vector = np.array([0., 0., -1.])

        # get in-plane z vector
        plane_normal = np.cross(in_plane_x_vector, in_plane_y_vector)
        plane_normal = plane_normal / np.linalg.norm(plane_normal)


        # resolve each contour in the list into the plane
        resolved_contours = []
        for contour in contours_list:
            # Translate the contour to the origin
            contour -= plane_origin
            # Rotate the contour to align with the in-plane x vector
            rotated_contour_x = np.dot(contour, in_plane_x_vector)
            rotated_contour_y = np.dot(contour, in_plane_y_vector)
            rotated_contour_z = np.dot(contour, plane_normal)

            # Create a new contour with the rotated coordinates
            resolved_contours.append(np.vstack([rotated_contour_x, rotated_contour_y, rotated_contour_z]).T)

        all_resolved_points = np.vstack(resolved_contours)
        all_resolved_points = np.unique(all_resolved_points, axis=0)

        # resolve the contour boundary into the plane
        resolved_boundary_x = np.dot(spline_contour_boundary, in_plane_x_vector)
        resolved_boundary_y = np.dot(spline_contour_boundary, in_plane_y_vector)
        resolved_boundary = np.vstack([resolved_boundary_x, resolved_boundary_y]).T

        resolved_spline_contours = []
        for contour in spline_contour_list:
            # Translate the contour to the origin
            contour -= plane_origin

            # Rotate the contour to align with the in-plane x vector
            rotated_contour_x = np.dot(contour, in_plane_x_vector)
            rotated_contour_y = np.dot(contour, in_plane_y_vector)

            # Create a new contour with the rotated coordinates
            resolved_spline_contours.append(np.vstack([rotated_contour_x, rotated_contour_y]).T)

            # refine end of spline contours by to match ends of fit contours
        
        refined_spline_contour_list = []
        for spline_contour, fit_contour in zip(resolved_spline_contours, resolved_contours):
            # sort by in_plane x coordinate
            spline_contour = spline_contour[np.argsort(spline_contour[:, 0])]   
            fit_contour = fit_contour[np.argsort(fit_contour[:, 0])]


            min_resolved_x_spline = min(spline_contour[:, 0])
            max_resolved_x_spline = max(spline_contour[:, 0])
            min_resolved_x_fit = min(fit_contour[:, 0])
            max_resolved_x_fit = max(fit_contour[:, 0])

            if min_resolved_x_spline < min_resolved_x_fit:
                # trim the spline contour
                spline_contour = spline_contour[spline_contour[:, 0] >= min_resolved_x_fit]
            if max_resolved_x_spline > max_resolved_x_fit:
                spline_contour = spline_contour[spline_contour[:, 0] <= max_resolved_x_fit]

            # Update the min and max resolved x values
            min_resolved_x_spline = min(spline_contour[:, 0])
            max_resolved_x_spline = max(spline_contour[:, 0])

            if min_resolved_x_spline > min_resolved_x_fit:
                new_start = fit_contour[fit_contour[:, 0] <= min_resolved_x_spline][:, :2]
                spline_contour = np.vstack([new_start, spline_contour])
            if max_resolved_x_spline < max_resolved_x_fit:
                new_end = fit_contour[fit_contour[:, 0] >= max_resolved_x_spline][:, :2]
                spline_contour = np.vstack([spline_contour, new_end])

            spline_contour = np.unique(spline_contour, axis=0)
            
            refined_spline_contour_list.append(spline_contour)

        all_resolved_spline_points = np.vstack(refined_spline_contour_list)

        interpolator = RBFInterpolator(all_resolved_points[:, :2], all_resolved_points[:, 2], kernel='thin_plate_spline')

        interpolated = interpolator(all_resolved_spline_points)

        reprojected_contours = []
        for contour in refined_spline_contour_list:
            interp_z = interpolator(contour)
            to_reproject = np.vstack([contour.T, interp_z]).T
            after_reprojection = plane_origin + to_reproject[:, 0][:, np.newaxis] * np.repeat([in_plane_x_vector], to_reproject.shape[0], axis=0) + to_reproject[:, 1][:, np.newaxis] * np.repeat([in_plane_y_vector], to_reproject.shape[0], axis=0) + to_reproject[:, 2][:, np.newaxis] * np.repeat([plane_normal], to_reproject.shape[0], axis=0)
            reprojected_contours.append(after_reprojection)
        reprojected_points = np.vstack(reprojected_contours)

        reprojected_contour_dict = {tuple(point): i for (i,point) in enumerate(reprojected_points)}

        reprojected_contour_boundary = np.vstack([reprojected_contours[0], 
                                    np.vstack([reprojected_contours[i][-1] for i in range(1, len(reprojected_contours) -1)]),
                                    reprojected_contours[-1][::-1],
                                    np.vstack([reprojected_contours[i][0] for i in range(len(reprojected_contours) -1, 0, -1)])])
        

                # resolve the contour boundary into the plane
        resolved_boundary_x = np.dot(reprojected_contour_boundary, in_plane_x_vector)
        resolved_boundary_y = np.dot(reprojected_contour_boundary, in_plane_y_vector)
        resolved_boundary = np.vstack([resolved_boundary_x, resolved_boundary_y]).T

        # Turn back into contour coordinates
        interpolated_xyz = plane_origin + all_resolved_spline_points[:, 0][:, np.newaxis] * np.repeat([in_plane_x_vector], interpolated.shape[0], axis=0) + all_resolved_spline_points[:, 1][:, np.newaxis] * np.repeat([in_plane_y_vector], interpolated.shape[0], axis=0) + interpolated[:, np.newaxis] * np.repeat([plane_normal], interpolated.shape[0], axis=0)
        if check_mesh:
            plt.close("all")
            fig, ax = plt.subplots()
            ax.plot(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", label="Resolved boundary")
            ax.scatter(resolved_boundary[:, 0], resolved_boundary[:, 1], color="blue", s=1, label="Resolved boundary points")
            plt.show(block=True)

        boundary_segments = np.array([[reprojected_contour_dict[tuple(reprojected_contour_boundary[i])], reprojected_contour_dict[tuple(reprojected_contour_boundary[i + 1])]] for i in range(len(reprojected_contour_boundary) - 1)])
        boundary_segments = np.vstack([boundary_segments, [reprojected_contour_dict[tuple(reprojected_contour_boundary[-1])], reprojected_contour_dict[tuple(reprojected_contour_boundary[0])]]])
        A = dict(vertices=all_resolved_spline_points, segments=boundary_segments)
        B = tr.triangulate(A, 'p')

        if check_mesh:
            plt.close("all")
            tr.compare(plt, A, B)
            plt.show(block=True)

        # write the mesh to a file
        mesh = meshio.Mesh(points=interpolated_xyz, cells={"triangle": B['triangles']})

        return mesh

    def mesh_simple_contours(self, depths: Union[np.ndarray, List[float]], mesh_name: str = None, mesh_format="vtk", check_mesh: bool=False, check_strike_dip: bool=True, km: bool = False):
        """
        Create a mesh from simple contours.

        :param depths: Array of depths for each contour
        :type depths: Union[np.ndarray, List[float]]
        :param km: Whether the depths are in kilometers, defaults to False
        :type km: bool
        :param mesh_name: Name of the output mesh file
        :type mesh_name: str, optional
        :param mesh_format: Format of the output mesh file
        :type mesh_format: str, optional
        :param check_mesh: Whether to display visualization of the mesh for debugging, defaults to False
        :type check_mesh: bool, optional
        :param check_strike_dip: Whether to print strike and dip information for debugging, defaults to False
        :type check_strike_dip: bool, optional
        :return: Triangular mesh representing the fault surface
        :rtype: meshio.Mesh
        """

        contours = self.generate_simple_contours(depths, km=km)

        # Create a mesh from the contours
        mesh = triangulate_contours(contours, mesh_name=mesh_name, mesh_format=mesh_format, check_mesh=check_mesh, check_strike_dip=check_strike_dip)

        return mesh
