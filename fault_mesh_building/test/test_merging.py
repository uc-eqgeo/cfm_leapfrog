from unittest import TestCase
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString

from fault_mesh.utilities.merging import align_two_nearly_adjacent_segments, merge_two_nearly_adjacent_segments, \
    merge_multiple_nearly_adjacent_segments, densify_line, sorted_merge

class test_merging(TestCase):
    def test_align_two_nearly_adjacent_segments(self):
        line1 = LineString([(0, 0), (1, 1)])
        line2 = LineString([(0, 1), (1, 0)])
        new_line1, new_line2 = align_two_nearly_adjacent_segments([line1, line2])
        self.assertTrue(new_line1.equals(LineString([(0, 0.5), (1, 1)])))
        self.assertTrue(new_line2.equals(LineString([(0, 0.5), (1, 0)])))

    def test_align_two_nearly_adjacent_segments_with_densify(self):
        pass


    def test_merge_two_nearly_adjacent_segments(self):
        line1 = LineString([(0, 0), (1, 1)])
        line2 = LineString([(0, 1), (1, 0)])
        new_line = merge_two_nearly_adjacent_segments([line1, line2])
        self.assertTrue(new_line.equals(LineString([(1, 0), (0, 0.5), (1, 1)])))

    def test_merge_two_nearly_adjacent_segments2(self):
        line1 = LineString([(0, 0), (200, 1)])
        line2 = LineString([(200, 0), (400, 0)])
        new_line = merge_two_nearly_adjacent_segments([line1, line2])
        self.assertTrue(new_line.equals(LineString([(0, 0), (200, 0.5), (400, 0)])))

    def test_merge_two_nearly_adjacent_segments_with_tolerance(self):
        line1 = LineString([(0, 0), (200, 1)])
        line2 = LineString([(200, 0), (400, 0)])
        new_line = merge_two_nearly_adjacent_segments([line1, line2], tolerance=1.)
        self.assertTrue(new_line.equals(LineString([(0, 0), (200, 0.5), (400, 0)])))

    def test_merge_two_nearly_adjacent_segments_with_tolerance2(self):
        line1 = LineString([(0, 0), (200, 1)])
        line2 = LineString([(200, 0), (400, 0)])
        with self.assertRaises(AssertionError):
            merge_two_nearly_adjacent_segments([line1, line2], tolerance=0.1)

    def test_merge_multiple_nearly_adjacent_segments(self):
        line1 = LineString([(0, 0), (1, 1)])
        line2 = LineString([(0, 1), (1, 0)])
        line3 = LineString([(0, 2), (1, 1)])
        new_line = merge_multiple_nearly_adjacent_segments([line1, line2, line3])
        self.assertTrue(new_line.equals(LineString([(0, 2), (1, 1), (0, 0.5), (1, 0)])))

    def test_densify_line(self):
        line = LineString([(0, 0), (1, 1)])
        new_line = densify_line(line, density=0.1)
        self.assertTrue(new_line.equals(LineString([(0, 0), (0.1, 0.1), (0.2, 0.2), (0.3, 0.3), (0.4, 0.4), (0.5, 0.5),
                                                    (0.6, 0.6), (0.7, 0.7), (0.8, 0.8), (0.9, 0.9), (1, 1)])))

    def test_sorted_merge(self):
        line1 = LineString([(0, 0), (1, 1)])
        line2 = LineString([(0, 1), (1, 0)])
        line3 = LineString([(0, 2), (1, 1)])
        new_line = sorted_merge([line1, line2, line3])
        self.assertTrue(new_line.equals(LineString([(0, 2), (1, 1), (0, 0.5), (1, 0)])))
