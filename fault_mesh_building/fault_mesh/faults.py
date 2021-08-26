from eq_fault_geom.geomio.cfm_faults import CfmMultiFault, CfmFault
import geopandas as gpd
import pandas as pd
from typing import List, Union
from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString

from fault_mesh.smoothing import straighten, smooth_trace


class LeapfrogMultiFault(CfmMultiFault):
    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, exclude_region_polygons: list = None,
                 exclude_region_min_sr: float = 1.8, include_names: list = None, depth_type: str = "D90",
                 exclude_aus: bool = True, exclude_zero: bool = True, sort_sr: bool = False):
        super(LeapfrogMultiFault, self).__init__(fault_geodataframe=fault_geodataframe, 
                                                 exclude_region_polygons=exclude_region_polygons,
                                                 exclude_region_min_sr=exclude_region_min_sr,
                                                 include_names=include_names,
                                                 depth_type=depth_type,
                                                 exclude_aus=exclude_aus, exclude_zero=exclude_zero, sort_sr=sort_sr)

    def add_fault(self, series: pd.Series, depth_type: str = "D90"):
        cfmFault = LeapfrogFault.from_series(series, parent_multifault=self, depth_type=depth_type)
        self.faults.append(cfmFault)


class LeapfrogFault(CfmFault):
    def __init__(self, parent_multifault: LeapfrogMultiFault = None, smoothing: int = None,
                 trimming_gradient: float = None):
        super(LeapfrogFault, self).__init__(None)
        self._parent = parent_multifault

        self.connections = []
        self._neighbouring_segments = None
        self._is_segment = False
        self._smoothing = smoothing
        self._trimming_gradient = trimming_gradient

    @property
    def is_segment(self):
        return self._is_segment

    @property
    def smoothing(self):
        return self._smoothing

    @property
    def trimming_gradient(self):
        return self._trimming_gradient

    @property
    def neighbouring_segments(self):
        return self._neighbouring_segments

    @neighbouring_segments.setter
    def neighbouring_segments(self, segment_list: list):
        assert isinstance(segment_list, list)
        assert len(segment_list) <= 2
        if len(segment_list) == 0:
            self._is_segment = False
        else:
            self._is_segment = True

        self._neighbouring_segments = segment_list

    def depth_contour(self, depth: float, smoothing: int = None, damping: int = None, km= False):
        if depth <= 0:
            shift = depth / self.down_dip_vector[-1]
        else:
            shift = (-1 * depth) / self.down_dip_vector[-1]
        if km:
            shift *= 1000.

        xo, yo, zo = shift * self.down_dip_vector
        contour = translate(self.nztm_trace, xo, yo, zo)

        if damping is not None:
            damped_contour = straighten(contour, self.dip_dir, damping)
        else:
            damped_contour = contour

        if smoothing is not None:
            smoothed_contour = smooth_trace(damped_contour, smoothing)

        elif self.smoothing is not None:
            smoothed_contour = smooth_trace(damped_contour, self.smoothing)

        else:
            smoothed_contour = damped_contour

        return smoothed_contour

    @property
    def nztm_trace(self):
        return self._nztm_trace

    @nztm_trace.setter
    def nztm_trace(self, trace: LineString):
        assert isinstance(trace, (LineString, MultiLineString))
        if isinstance(trace, MultiLineString):
            trace = list(trace)[0]

        new_trace = LineString([(xi, yi, 0.) for xi, yi in trace.coords])

        self._nztm_trace = new_trace







class LeapfrogConnectedFault:
    def __init__(self, fault_list: List[LeapfrogFault], parent_multifault: LeapfrogMultiFault = None):
        assert all([isinstance(x, LeapfrogFault) for x in fault_list])
        self._parent = parent_multifault
        self.faults = fault_list

        self.connections = []




    def join_combine_surface_traces(self, smoothing: int = 0, simplification: int = 0):
        pass

    def depth_contours(self, smoothing: int = 0, simplification: int = 0, angle_gap: Union[float, int] = None):
        pass


class LeapfrogFaultSegment(LeapfrogFault):
    def __init__(self, parent_fault: LeapfrogConnectedFault, parent_multifault: LeapfrogMultiFault = None, ):
        super(LeapfrogFault, self).__init__(parent_multifault)
        self._fault = parent_fault
        self.neighbouring_segments = []




class LeapfrogFaultModel:
    def __init__(self, fault_list: List[Union[LeapfrogFault, LeapfrogConnectedFault]]):
        assert all([isinstance(x, (LeapfrogFault, LeapfrogConnectedFault)) for x in fault_list])
        self.faults = fault_list
