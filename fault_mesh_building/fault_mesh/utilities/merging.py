"""
Functions for merging segments that are nearly adjacent.
"""

import numpy as np
from typing import List
from shapely.geometry import LineString, Point
from shapely.ops import linemerge


def align_two_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200.):
    assert len(segment_list) == 2
    line1, line2 = segment_list
    assert line1.distance(line2) <= tolerance

    l1e1 = Point(line1.coords[0])
    l1e2 = Point(line1.coords[-1])
    p1 = l1e1 if l1e1.distance(line2) <= l1e2.distance(line2) else l1e2

    l2e1 = Point(line2.coords[0])
    l2e2 = Point(line2.coords[-1])

    p2 = l2e1 if l2e1.distance(line1) <= l2e2.distance(line1) else l2e2

    mid_point = Point(0.5 * (np.array(p1.coords) + np.array(p2.coords)).flatten())

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
        while len(segment_list) > 2:
            segment_list = [merge_two_nearly_adjacent_segments(segment_list[:2], tolerance)] + segment_list[2:]
        return merge_two_nearly_adjacent_segments(segment_list, tolerance)