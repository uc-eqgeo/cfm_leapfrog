"""
Functions for merging segments that are nearly adjacent.
This module provides functions to align and merge line segments that are close to each other.
It includes functions to densify line segments, align two nearly adjacent segments, merge two nearly adjacent segments, and merge multiple nearly adjacent segments. The merging process is based on the distance between the segments and can be customized with a tolerance parameter.
The functions are designed to work with Shapely LineString objects and can be used in geospatial applications where line segments need to be combined or aligned.   

Functions:
    - densify_line: Densifies a LineString by adding interpolated points along the line.
    - align_two_nearly_adjacent_segments: Aligns two nearly adjacent line segments by moving their closest endpoints to a common midpoint.
    - merge_two_nearly_adjacent_segments: Merges two nearly adjacent line segments into a single line.
    - merge_multiple_nearly_adjacent_segments: Merges multiple line segments that are nearly adjacent to each other.
    - sorted_merge: Merges a list of LineString segments into a single LineString by iteratively combining the nearest segments.
"""

import numpy as np
from typing import List, Union
from shapely.geometry import LineString, Point
from shapely.ops import linemerge

def densify_line(line: LineString, density: float = 1.3):
    """
    Densify a LineString by adding interpolated points along the line.
    This function takes a LineString and adds points at regular intervals
    determined by the density parameter. The original points are preserved,
    and the result is a LineString with a higher density of points.
    :param line: The input LineString to be densified
    :type line: LineString
    :param density: The interval between interpolated points along the line, defaults to 1.3
    :type density: float, optional
    :return: A new LineString with additional interpolated points
    :rtype: LineString
    """

    original_coords = [Point(*coord) for coord in line.coords]
    interp_dists = np.arange(density, line.length, density)
    interpolated_points = [line.interpolate(dist) for dist in interp_dists]
    combined_coords = original_coords + interpolated_points
    sorted_coords = sorted(combined_coords, key=lambda x: line.project(x))
    return LineString(sorted_coords)

def align_two_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200., densify: float = 1.e3):
    """Aligns two nearly adjacent line segments by moving their closest endpoints to a common midpoint.
    This function takes two LineString objects that are within a specified distance tolerance of each other
    and aligns them by replacing their closest endpoints with a common midpoint. The function can optionally
    densify the line segments before alignment to improve the precision of the closest-point determination.
    :param segment_list: List containing exactly two LineString objects that need to be aligned
    :param tolerance: Maximum distance between segments to be considered for alignment, defaults to 200.
    :param densify: Distance between points when densifying the lines, defaults to 1.e3.
                   If None, lines are not densified.
    :return: A tuple containing two aligned LineString objects
    :rtype: tuple(LineString, LineString)
    :raises AssertionError: If segment_list doesn't contain exactly 2 LineStrings
    :raises AssertionError: If the distance between the two LineStrings exceeds the tolerance
    """
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

def merge_two_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200., densify: float = 1.e3):
    """Merge two nearly adjacent line segments into a single line.
    This function takes two line segments that are close to each other (within a specified tolerance)
    and merges them into a single continuous LineString. It first aligns the segments using
    the align_two_nearly_adjacent_segments function, and then merges them.
    :param segment_list: List containing two LineString objects to be merged
    :param tolerance: Maximum distance between endpoints for segments to be considered adjacent, defaults to 200.
    :param densify: Distance between points when densifying lines for better alignment, defaults to 1.e3
    :return: A single merged LineString created from the two input segments
    :rtype: LineString
    """
    new_line1, new_line2 = align_two_nearly_adjacent_segments(segment_list, tolerance, densify=densify)
    return linemerge([new_line1, new_line2])


def merge_multiple_nearly_adjacent_segments(segment_list: List[LineString], tolerance: float = 200., densify: Union[float, None] = 1.e3):
    """Merge multiple line segments that are nearly adjacent to each other.
    This function takes a list of LineString objects and merges them into a single LineString
    if they are within the specified tolerance distance of each other. If there are only
    two segments, it uses merge_two_nearly_adjacent_segments; otherwise, it uses sorted_merge.
    :param segment_list: List of LineString objects to be merged
    :param tolerance: Maximum distance between segments to be considered for merging, defaults to 200.0
    :param densify: Distance between points when densifying the lines, set to None to skip densification, defaults to 1.0e3
    :return: A merged LineString containing all input segments
    :rtype: LineString                                  
    """
    assert len(segment_list) >= 2
    if len(segment_list) == 2:
        return merge_two_nearly_adjacent_segments(segment_list, tolerance, densify=densify)
    else:

        return sorted_merge(segment_list, tolerance, densify=densify)

def sorted_merge(seg_list: List[LineString], tolerance: float = 200., densify: float = 1.e3):
    """Merges a list of LineString segments into a single LineString by iteratively combining the nearest segments.
    This function starts with the southernmost segment (minimum y-coordinate of centroid) and 
    repeatedly merges it with the nearest remaining segment until all segments are included.
    :param seg_list: List of LineString segments to be merged
    :param tolerance: Maximum distance allowed between segments for merging, defaults to 200.
    :param densify: Distance between points when densifying segments for better merging results, defaults to 1.e3
    :return: A single LineString representing the merged segments
    :rtype: LineString
    """
    sorted_list = [min(seg_list, key=lambda x:x.centroid.y)]
    combined_segment = sorted_list[0]
    while len(sorted_list) < len(seg_list):
        remaining_list = [seg for seg in seg_list if seg not in sorted_list]
        nearest_remaining = sorted(remaining_list, key=lambda x: x.distance(combined_segment))[0]
        combined_segment = merge_two_nearly_adjacent_segments([combined_segment, nearest_remaining],
                                                              tolerance=tolerance, densify=densify)
        sorted_list.append(nearest_remaining)

    return combined_segment