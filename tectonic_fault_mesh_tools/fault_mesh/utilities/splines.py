"""Utilities for working with splines and contour interpolation.
This module provides functions for fitting splines to 2D contour data,
creating linearly spaced arrays, and smoothing contour geometries.
It is particularly useful for processing geological fault data and
other linear features that require smoothing or regularization.
Functions:
    fit_2d_line: Calculate dip angle by fitting a line to 2D data points
    linspace_with_spacing: Create evenly spaced points with a specified spacing
    spline_fit_contours: Fit smooth splines to contour geometries in a GeoDataFrame
"""
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from scipy.interpolate import make_interp_spline, UnivariateSpline

def fit_2d_line(x: np.ndarray, y: np.ndarray):
    """
    Fit a 2D line to a set of points
    :param x:
    :param y:
    :return:
    """
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)

    px = np.polyfit(x, y, 1, full=True)
    gradient_x = px[0][0]

    if len(px[1]):
        res_x = px[1][0]
    else:
        res_x = 0

    py = np.polyfit(y, x, 1, full=True)
    gradient_y = py[0][0]
    if len(py[1]):
        res_y = py[1][0]
    else:
        res_y = 0

    if res_x <= res_y:
        dip_angle = np.degrees(np.arctan(gradient_x))
    else:
        dip_angle = np.degrees(np.arctan(1./gradient_y))

    return dip_angle


def linspace_with_spacing(start: float, stop: float, spacing: float) -> np.ndarray:
    """
    Create a linearly spaced array with a given spacing
    :param start:
    :param stop:
    :param spacing:
    :return:
    """
    num_points = int(np.ceil((stop - start) / spacing)) + 1
    return np.linspace(start, stop, num_points, endpoint=True)

def spline_fit_contours(contours: gpd.GeoDataFrame, point_spacing: float = 100., output_spacing: float = 1000.) -> gpd.GeoDataFrame:
    """Fits a smooth spline to each contour in a GeoDataFrame.
    This function takes contour geometries, fits a spline through them, and returns a new set of smoothed contours
    with points at a regular spacing. The process involves segmenting each contour, transforming it to align with 
    its principal direction (strike), fitting a spline in this transformed space, and then transforming back to 
    the original coordinate system.
    :param contours: GeoDataFrame containing LineString geometries representing the contours to be smoothed.
    :param point_spacing: Distance between points when segmenting the original contour for spline fitting, defaults to 100.
    :param output_spacing: Distance between points in the resulting smoothed contour, defaults to 1000.
    :return: GeoDataFrame containing the smoothed contours as LineString geometries.
    """
    # Create a list to store interpolated contours
    interpolated_contours = []

    # Iterate over each contour
    for contour in contours.geometry:
        # Segmentize the contour
        segmentized_contour = contour.segmentize(point_spacing)
        # Convert the segmentized contour to a numpy array
        try:
            segmentized_contour = np.vstack([np.array(segment.coords) for segment in segmentized_contour.geoms])
        except AttributeError:
            # If the contour is not a MultiLineString, convert it directly
            segmentized_contour = np.array(segmentized_contour.coords)
        # Find the overall strike of the contour
        strike = 90 - fit_2d_line(segmentized_contour[:, 0], segmentized_contour[:, 1])
        # Create a strike vector
        strike_vector = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.])
        strike_vector /= np.linalg.norm(strike_vector)
        # Across strike vector
        across_strike_vector = np.array([-strike_vector[1], strike_vector[0], 0.])
        # centroid of the contour
        centroid = np.mean(segmentized_contour, axis=0)
        # Translate the contour to the origin
        segmentized_contour -= centroid
        # rotate the contour to align with the strike vector
        rotated_contour_x = np.dot(segmentized_contour, strike_vector)
        rotated_contour_y = np.dot(segmentized_contour, across_strike_vector)
        # Create a new contour with the rotated coordinates
        rotated_contour = np.vstack([rotated_contour_x, rotated_contour_y]).T
        # sort the contour by the rotated x coordinate
        rotated_contour = rotated_contour[np.argsort(rotated_contour[:, 0])]

        # Remove duplicate points before spline fitting
        unique_mask = np.concatenate(([True], np.any(np.diff(rotated_contour, axis=0) != 0, axis=1)))
        rotated_contour = rotated_contour[unique_mask]

        # Create a spline for the rotated contour
        spline = make_interp_spline(rotated_contour[:, 0], rotated_contour[:, 1], k=2)
        # Append the spline to the list
        interp_x = linspace_with_spacing(rotated_contour[0, 0], rotated_contour[-1, 0], output_spacing)
        interp_y = np.array([spline(x) for x in interp_x])
        # Create a new contour with the interpolated coordinates
        interp_contour = np.vstack([interp_x, interp_y]).T
        # Rotate the contour back to the original coordinates
        rotated_contour = np.column_stack([interp_contour[:, 0] * strike_vector[0] + interp_contour[:, 1] * across_strike_vector[0],
                                          interp_contour[:, 0] * strike_vector[1] + interp_contour[:, 1] * across_strike_vector[1],
                                          np.zeros_like(interp_contour[:, 0])])
        # Translate the contour back to the original coordinates
        rotated_contour += centroid
        # Append the rotated contour to the list of interpolated contours
        interpolated_contours.append(LineString(rotated_contour))
        # Create a new GeoDataFrame with the interpolated contours
        interpolated_contours_gdf = gpd.GeoDataFrame(geometry=interpolated_contours, crs=contours.crs)
    return interpolated_contours_gdf


def fault_edge_spline(max_slip: float, distance_width: float, total_width: float, min_slip_fraction: float = 0., line_stop_fraction: float = 0.5, gradient_change: float = 1.2,
                      spline_k: int = 3, spline_s: float = 0.0, resolution: float = 10.):
    """
    Create a spline that can be used to create a slip distribution along the edge of a fault.
    @param max_slip:
    @param distance_width:
    @param min_slip_fraction:
    @param line_stop_fraction:
    @param gradient_change:
    @param spline_k:
    @param spline_s:
    @param resolution:
    @return:
    """

    x_points = np.hstack([np.arange(0., distance_width * line_stop_fraction + resolution, resolution),
                            np.arange(distance_width, total_width + resolution, resolution)])

    x_for_interp = np.array([0., distance_width * line_stop_fraction, distance_width, total_width])
    y_for_interp = np.array([max_slip * min_slip_fraction, gradient_change *
                             (max_slip * min_slip_fraction + line_stop_fraction * max_slip * (1 - min_slip_fraction)),
                             max_slip, max_slip])

    interpolated_y = np.interp(x_points, x_for_interp, y_for_interp)
    out_spline = UnivariateSpline(x_points, interpolated_y, k=spline_k, s=spline_s)
    return out_spline

def fault_depth_spline(gradient_change_x: float, after_change_fract: float, resolution: float = 0.01,
                       after_change_gradient: float = 1.2, spline_k: int = 3, spline_s: float = 0.0):
    x_points = np.hstack([np.arange(0., gradient_change_x + resolution, resolution),
                            np.arange(1 - (1. - gradient_change_x) * after_change_fract, 1. + resolution, resolution)])
    x_for_interp = np.array([0., gradient_change_x, 1. - (1. - gradient_change_x) * after_change_fract, 1.])
    y_for_interp = np.array([1., 1., after_change_gradient * after_change_fract, 0.])
    interpolated_y = np.interp(x_points, x_for_interp, y_for_interp)
    out_spline = UnivariateSpline(x_points, interpolated_y, k=spline_k, s=spline_s)
    return out_spline

