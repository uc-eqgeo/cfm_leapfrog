from typing import List
from shapely.geometry import LineString, Point


def cut(line: LineString, distance: float):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line), LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i + 1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            if line.has_z:
                return LineString(coords[:i] + [(cp.x, cp.y, cp.z)]), LineString([(cp.x, cp.y, cp.z)] + coords[i:])
            else:
                return LineString(coords[:i] + [(cp.x, cp.y)]), LineString([(cp.x, cp.y)] + coords[i:])


def cut_line_at_point(line: LineString, point: Point):
    assert isinstance(point, Point)
    distance = line.project(point)
    l1, l2 = cut(line, distance)
    return l1, l2


def cut_line_at_multiple_points(line: LineString, points: List[Point]):
    sorted_points = sorted(points, key=lambda x: line.project(x))
    segments = []
    line_to_cut = None
    for point in sorted_points[:-1]:
        if line_to_cut is None:
            l1, l2 = cut_line_at_point(line, point)
        else:
            l1, l2 = cut_line_at_point(line_to_cut, point)
        line_to_cut = l2
        segments.append(l1)
    l1, l2 = cut_line_at_point(line_to_cut, sorted_points[-1])
    segments += [l1, l2]

    return segments


def cut_line_between_two_points(line: LineString, points: List[Point]):
    assert len(points) == 2
    segments = cut_line_at_multiple_points(line, points)
    return segments[1]
