from typing import Union
import fnmatch
from shapely.geometry import MultiLineString
from shapely.ops import linemerge

from fault_mesh.smoothing import smooth_trace
from fault_mesh.utilities.cutting import cut_line_at_multiple_points

class ConnectedFaultSystem:
    def __init__(self, overall_name: str, search_patterns: Union[str, list], cfm_faults,
                 excluded_names: Union[str, list] = None, tolerance: float = 100., smooth_trace_refinements: int = None):
        self.overall_name = overall_name
        segments = []

        if isinstance(search_patterns, str):
            search_pattern_list = [search_patterns]
        else:
            assert isinstance(search_patterns, list)
            assert len(search_patterns)
            search_pattern_list = search_patterns

        if isinstance(excluded_names, str):
            excluded_list = [excluded_names]
        elif excluded_names == None:
            excluded_list = []
        else:
            assert isinstance(excluded_names, list)
            assert len(excluded_names)
            excluded_list = excluded_names

        for fault in cfm_faults.faults:
            name = fault.name
            if any([fnmatch.fnmatch(name, pattern) for pattern in search_pattern_list]):
                if not any([fnmatch.fnmatch(name, pattern) for pattern in excluded_list]):
                    segments.append(fault)

        for segment in segments:
            other_segments = MultiLineString([other_seg.nztm_trace for other_seg in segments if other_seg != segment])
            closest_dist = segment.nztm_trace.distance(other_segments)
            if closest_dist > tolerance:
                raise ValueError(f"Fault traces >{tolerance} m apart: {segment.name}")

        # Order segments
        # Find furthest south or west

        sorted_segments_s_to_n = sorted(segments, key = lambda x: (x.nztm_trace.bounds[0],
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

        self.segment_boundaries = [l1.nztm_trace.intersection(l2.nztm_trace) for l1, l2 in zip(self.segments,
                                                                                               self.segments[1:])]

        self.overall_trace = linemerge([seg.nztm_trace for seg in self.segments])

        if smooth_trace_refinements is not None:
            self.smoothed_overall_trace = smooth_trace(self.overall_trace, n_refinements=smooth_trace_refinements)
            self.smoothed_segments = cut_line_at_multiple_points(self.smoothed_overall_trace, self.segment_boundaries)
        else:
            self.smoothed_overall_trace = None
            self.smoothed_segments = None





        # Might need to add interpolation between nearby but not identical segment ends. Do this when we get to it.










