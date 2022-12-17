"""
Functions for merging segments that are nearly adjacent.
"""

import numpy as np
from typing import List
from shapely.geometry import LineString, Point
from shapely.ops import linemerge

def densify_line(line: LineString, density: float = 1.3):
    original_coords = [Point(*coord) for coord in line.coords]
    interp_dists = np.arange(density, line.length, density)
    interpolated_points = [line.interpolate(dist) for dist in interp_dists]
    combined_coords = original_coords + interpolated_points
    sorted_coords = sorted(combined_coords, key=lambda x: line.project(x))
    return LineString(sorted_coords)

def align_two_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200., densify: float = 1.e3):
    assert len(segment_list) == 2
    line1, line2 = segment_list
    assert line1.distance(line2) <= tolerance

    if densify is not None:
        line1 = densify_line(line1, density=densify)
        line2 = densify_line(line2, density=densify)

    l1e1 = Point(line1.coords[0])
    l1e2 = Point(line1.coords[-1])
    p1 = l1e1 if l1e1.distance(line2) <= l1e2.distance(line2) else l1e2

    l2e1 = Point(line2.coords[0])
    l2e2 = Point(line2.coords[-1])

    p2 = l2e1 if l2e1.distance(line1) <= l2e2.distance(line1) else l2e2

    mid_point = Point(0.5 * (np.array(p1.coords[0]) + np.array(p2.coords[0])))

    if l1e1 == p1:
        new_line1 = np.vstack([np.array(mid_point.coords), np.array(line1.coords)[1:]])
    else:
        new_line1 = np.vstack([np.array(line1.coords)[:-1], np.array(mid_point.coords)])

    if l2e1 == p2:
        new_line2 = np.vstack([np.array(mid_point.coords), np.array(line2.coords)[1:]])
    else:
        new_line2 = np.vstack([np.array(line2.coords)[:-1], np.array(mid_point.coords)])

    return LineString(new_line1), LineString(new_line2)

def merge_two_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200.):
    new_line1, new_line2 = align_two_nearly_adjacent_segments(segment_list, tolerance)
    return linemerge([new_line1, new_line2])


def merge_multiple_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200.):
    assert len(segment_list) >= 2
    if len(segment_list) == 2:
        return merge_two_nearly_adjacent_segments(segment_list, tolerance)
    else:

        return sorted_merge(segment_list, tolerance)

def sorted_merge(seg_list: List[LineString], tolerance: float = 200.):
    sorted_list = [min(seg_list, key=lambda x:x.centroid.y)]
    combined_segment = sorted_list[0]
    while len(sorted_list) < len(seg_list):
        remaining_list = [seg for seg in seg_list if seg not in sorted_list]
        nearest_remaining = sorted(remaining_list, key=lambda x: x.distance(combined_segment))[0]
        combined_segment = merge_two_nearly_adjacent_segments([combined_segment, nearest_remaining],
                                                              tolerance=tolerance)
        sorted_list.append(nearest_remaining)

    return combined_segment