from fault_mesh.faults.generic import (smallest_difference, normalize_bearing, bearing_leq,
                                       bearing_geq, reverse_bearing, reverse_line, calculate_dip_direction)

from unittest import TestCase
from shapely.geometry import LineString
import numpy as np


class test_bearing_functions(TestCase):
    def test_smallest_difference1(self):
        small_diff = smallest_difference(5, 355)
        self.assertAlmostEqual(small_diff, 10.)

    def test_smallest_difference2(self):
        small_diff = smallest_difference(-4., 355)
        self.assertAlmostEqual(small_diff, 1.)

    def test_smallest_difference3(self):
        small_diff = smallest_difference(46., 225)
        self.assertAlmostEqual(small_diff, 179)

    def test_smallest_difference4(self):
        small_diff = smallest_difference(44., 225)
        self.assertAlmostEqual(small_diff, 179)

    def test_normalize1(self):
        self.assertAlmostEqual(normalize_bearing(-45), 315)

    def test_normalize2(self):
        self.assertAlmostEqual(normalize_bearing(693), 333)

    def test_leq1(self):
        self.assertIs(bearing_leq(271., 90.), True)

    def test_leq2(self):
        self.assertIs(bearing_leq(269., 90.), False)

    def test_geq1(self):
        self.assertIs(bearing_geq(140., 315.), False)

    def test_geq2(self):
        self.assertIs(bearing_geq(90., 315.), True)

    def test_reverse_bearing1(self):
        self.assertAlmostEqual(180.1, reverse_bearing(0.1))



class test_geometric_functions(TestCase):
    def setUp(self) -> None:
        self.straight_fault_array = np.array([[1640000., 5350000.], [1650000., 5360000.]])
        self.array_reversed = self.straight_fault_array[-1::-1, :]
        self.straight_fault_linestring = LineString(self.straight_fault_array)
        self.bent_fault_array = np.array([[1640000., 5350000.], [1650000., 5360000.],
                                          [1640000., 5370000.]])
        self.bent_array_reversed = self.bent_fault_array[-1::-1, :]
        self.bent_line = LineString(self.bent_array_reversed)

    def test_reverse_line1(self):
        reversed_straight = reverse_line(self.straight_fault_linestring)
        np.testing.assert_array_almost_equal(np.array(reversed_straight.coords), self.array_reversed)

    def test_reverse_line2(self):
        reversed_bent = reverse_line(self.bent_line)
        np.testing.assert_array_almost_equal(np.array(reversed_bent.coords), self.bent_fault_array)

    def test_dip_direction_straight1(self):
        self.assertAlmostEqual(calculate_dip_direction(self.straight_fault_linestring), 135.)

    def test_dip_direction_straight2(self):
        self.assertAlmostEqual(calculate_dip_direction(LineString(reverse_line(self.straight_fault_linestring))), 315)

    def test_dip_direction_bent1(self):
        self.assertAlmostEqual(calculate_dip_direction(self.bent_line), 270.)

    def test_dip_direction_bent2(self):
        self.assertAlmostEqual(calculate_dip_direction(LineString(self.bent_fault_array)), 90.)

# class test_dip_direction(TestCase):








