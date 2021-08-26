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


def straighten(line: LineString, strike: float, damping: float):
    strike_vector = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.])
    across_strike = np.array([np.sin(np.radians(strike + 90.)), np.cos(np.radians(strike + 90.)), 0.])
    line_array = np.array(line)
    centroid = np.array(line.centroid)

    along_dists = np.dot(line_array - centroid, strike_vector)
    across_dists = np.dot(line_array - centroid, across_strike)

    new_locations = centroid + along_dists + damping * across_dists

    return LineString(new_locations)



