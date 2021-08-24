import numpy as np
from typing import Union
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString



def chaikins_corner_cutting(coords, refinements=5):
    coords = np.array(coords)

    for _ in range(refinements):
        l = coords.repeat(2, axis=0)
        R = np.empty_like(l)
        R[0] = l[0]
        R[2::2] = l[1:-1:2]
        R[1:-1:2] = l[2::2]
        R[-1] = l[-1]
        coords = l * 0.75 + R * 0.25

    return coords

def smooth_trace(trace: LineString, n_refinements: int = 5):
    assert isinstance(trace, LineString)
    coords = np.array(trace.coords)
    return LineString(chaikins_corner_cutting(coords, refinements=n_refinements))
