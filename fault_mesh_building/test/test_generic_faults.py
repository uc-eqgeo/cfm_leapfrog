from unittest import TestCase
import geopandas as gpd
import logging
import numpy as np
from numpy.testing import assert_array_almost_equal

from fault_mesh.faults.generic import GenericFault, GenericMultiFault, required_fields, expected_fields
import os

file_path = os.path.dirname(os.path.realpath(__file__))




class test_generic_faults(TestCase):

    def setUp(self):

        self.filename = os.path.join(file_path, "./data/selected_faults.shp")
        self.fault_geodataframe = gpd.read_file(self.filename)
        self.fault_model = GenericMultiFault(self.fault_geodataframe)
        self.logger = logging.getLogger('fault_model_logger')

        # Sort alphabetically by name
        self.sorted_df = self.fault_geodataframe.sort_values("Name")
        # Reset index to line up with alphabetical sorting
        self.sorted_df = self.sorted_df.reset_index(drop=True)

        self.faults = []


    def test_check_input1(self):
        df_response = self.fault_geodataframe[required_fields[:-1]].copy()
        with self.assertRaises(ValueError):
            self.fault_model.check_input1(df_response)

    def test_check_input2(self):
        exf = [i for i in expected_fields if i not in ['Depth_min', 'Name']]
        df_response = self.fault_geodataframe[exf[:-1]].copy()
        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.fault_model.check_input2(df_response)
            print(cm.output)
            self.assertIn(
                "WARNING:fault_model_logger:missing expected field", cm.output[0]
            )

    def test_add_fault(self):
        self.assertGreater(len(self.sorted_df), 0, "Not enough rows in the input data set")

        for i, fault in self.sorted_df.iterrows():
            length = len(self.fault_model.faults)
            self.fault_model.add_fault(fault)
            self.assertAlmostEqual(length+1, len(self.fault_model.faults))


    def test_fault_numbers(self):
        if self.assertIsNotNone(self.fault_model.fault_numbers):
            self.assertIsInstance(self.fault_model.fault_numbers, int)

        self.assertFalse(len(self.fault_model.fault_numbers) == 0, 'The fault number is missing')

    #assert False

    def test_from_shp(self):
        multi_fault = self.fault_model.from_nz_cfm_shp(self.filename)
        self.assertIsNotNone(multi_fault)
        response = isinstance(multi_fault, GenericMultiFault)
        self.assertTrue(response, 'supplied object is not a GenericMultiFault"'
                                  ', it is a "{}"'.format(type(multi_fault)))




class test_generic_fault(TestCase):
    def setUp(self):
        self.cmf_fault = GenericFault()
        self.logger = logging.getLogger('fault_model_logger')

        self.filename = os.path.join(file_path, "./data/selected_faults.shp")
        self.fault_geodataframe = gpd.GeoDataFrame.from_file(self.filename)
        self.fault_model = GenericMultiFault(self.fault_geodataframe)
        # Sort alphabetically by name
        self.sorted_df = self.fault_geodataframe.sort_values("Name")
        # Reset index to line up with alphabetical sorting
        self.sorted_df = self.sorted_df.reset_index(drop=True)




    def test_dip_max(self):
        self.cmf_fault.dip_max = 10.5
        self.assertAlmostEqual(self.cmf_fault.dip_max, 10.5)

        self.cmf_fault._dip_max = 8.6
        self.assertNotEqual(self.cmf_fault.dip_max, 10.5)
        #
        with self.assertRaises(Exception):
            self.cmf_fault.dip_max = "Hello"


        self.cmf_fault.dip_min = 20.6
        self.cmf_fault.dip_best = 40.1
        dip = 19.5

        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.cmf_fault.dip_max = dip
            self.assertIn(
                "WARNING:fault_model_logger:dip_max is lower than dip min or dip best", cm.output
            )


    def test_dip_min(self):
        self.cmf_fault.dip_min = 10.5
        self.assertAlmostEqual(self.cmf_fault.dip_min, 10.5)

        self.cmf_fault._dip_min = 8.6
        self.assertNotEqual(self.cmf_fault.dip_min, 10.5)
        #
        with self.assertRaises(Exception):
            self.cmf_fault.dip_min = "Hello"

        self.cmf_fault.dip_max = 45.3
        self.cmf_fault.dip_best = 40.1
        dip = 50.6

        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.cmf_fault.dip_min = dip
            self.assertIn(
                "WARNING:fault_model_logger:dip_min is higher than dip max or dip best", cm.output
            )



#not sure if the test beolw is correct;
    def test_dip_dir_str(self):
        dip_dir = 'SW'
        self.cmf_fault.dip_dir_str = dip_dir
        self.assertIsInstance(dip_dir, str)

        series = self.sorted_df.iloc[1]
        self.cmf_fault.nztm_trace = series['geometry']

        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.cmf_fault.validate_dip_direction()
            self.assertIn(
                "WARNING:fault_model_logger:Supplied trace and dip direction are inconsistent", cm.output
            )

        dip_dir = None
        self.cmf_fault.dip_dir_str = dip_dir
        self.assertAlmostEqual(self.cmf_fault.dip_dir, 337.261792201337)


    #still working on this
    def test_dip_sigma(self):
        self.cmf_fault._dip_sigma = 8.6
        self.assertAlmostEqual(self.cmf_fault.dip_sigma, 8.6)

        self.cmf_fault._dip_sigma = None
        self.cmf_fault.dip_min = 16
        self.cmf_fault.dip_max = 25
        self.assertAlmostEqual(self.cmf_fault.dip_sigma, 4.5)

    def test_validate_dip_direction(self):
        series = self.sorted_df.iloc[0]
        self.cmf_fault.nztm_trace = series['geometry']

        dip_dir = 'SE'
        self.cmf_fault.dip_dir_str = dip_dir

        self.cmf_fault.validate_dip_direction()
        self.assertAlmostEqual(self.cmf_fault.dip_dir, 165.33461826424696)

        dip_dir = None
        self.cmf_fault.dip_dir_str = dip_dir
        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.cmf_fault.validate_dip_direction()
            self.assertIn(
                "WARNING:fault_model_logger:Insufficient information to validate dip direction", cm.output
            )

        dip_dir = 'NE'
        self.cmf_fault.dip_dir_str = dip_dir
        with self.assertLogs(logger=self.logger, level='WARNING') as cm:
            self.cmf_fault.validate_dip_direction()
            self.assertIn(
                "WARNING:fault_model_logger:Supplied trace and dip direction are inconsistent", cm.output
            )


    def test_validate_dip(self):
        dip = 15.6
        self.assertIsInstance(self.cmf_fault.validate_dip(dip), float)

        dip = "Hello"
        with self.assertRaises(Exception):
            self.cmf_fault.validate_dip(dip)

        dip = -20.6  # should be between 0 - 90 otherwise assert error
        with self.assertRaises(Exception):
            self.cmf_fault.validate_dip(dip)

    def test_nztm_trace(self):
        series = self.sorted_df.iloc[0]
        trace = series['geometry']
        #trace = 0.124
        self.cmf_fault.nztm_trace = trace

    def test_wgs_trace(self):
        series = self.sorted_df.iloc[0]
        trace = series['geometry']
        self.cmf_fault.nztm_trace = trace

        reponseX, reponseY = self.cmf_fault.wgs_trace.coords.xy
        response = np.array(reponseX.tolist())
        actual = np.array([173.81651055, 173.75, 173.6, 173.517, 173.428])
        assert_array_almost_equal(response, actual)




