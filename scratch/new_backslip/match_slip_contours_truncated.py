import geopandas as gpd
from shapely.geometry import LineString
import numpy as np
from scipy.interpolate import make_interp_spline, UnivariateSpline

def spline_points_start(x_start: float,transition_dist: float, y_scaling: float, line_stop_frac: float = 0.5, before_gradient = 1.2, resolution = 200.):
    """
    Generate spline points for the start of a fault segment.
    
    Parameters:
    - x_start: Starting x-coordinate.
    - transition_dist: Distance for the transition.
    - y_scaling: Scaling factor for the y-coordinate.
    - line_stop_frac: Fraction of the line to stop at.
    - before_gradient: Gradient before the transition.
    - resolution: Resolution for the spline points.
    
    Returns:
    - Spline points as a tuple of x and y coordinates.
    """
    x_points = np.arange(x_start, x_start + transition_dist * line_stop_frac + resolution, resolution)
    y_points = np.interp(x_points, [x_points[0], x_points[-1]], [0, before_gradient * line_stop_frac])
    y_points *= y_scaling

    return x_points, y_points

def spline_points_end(x_end: float, transition_dist: float, y_scaling: float, line_stop_frac: float = 0.5, after_gradient = 1.2, resolution = 200.):
    """
    Generate spline points for the end of a fault segment.
    
    Parameters:
    - x_end: Starting x-coordinate.
    - transition_dist: Distance for the transition.
    - y_scaling: Scaling factor for the y-coordinate.
    - line_stop_frac: Fraction of the line to stop at.
    - after_gradient: Gradient after the transition.
    - resolution: Resolution for the spline points.
    
    Returns:
    - Spline points as a tuple of x and y coordinates.
    """
    x_points = np.arange(x_end - transition_dist * line_stop_frac, x_end + resolution, resolution)
    y_points = y_scaling * np.interp(x_points, [x_points[0], x_points[-1]], [after_gradient * line_stop_frac, 0])

    return x_points, y_points

def spline_points_middle(x_start: float, x_end: float, transition_dist1: float, transition_dist2: float, 
                         y_start: float, y_end: float, resolution: float = 200., gradient: float = 1.2,
                         line_stop_frac: float = 0.5):
    """
    Generate spline points for the middle of a fault segment.

    Parameters:
    - x_start: Starting x-coordinate.
    - x_end: Ending x-coordinate.
    - transition_dist1: Distance for the first transition.
    - transition_dist2: Distance for the second transition.
    - y_start: Starting y-coordinate.
    - y_end: Ending y-coordinate.
    - resolution: Resolution for the spline points.
    - gradient: Gradient for the spline points.
    - line_stop_frac: Fraction of the line to stop at.

    Returns:
    - Spline points as a tuple of x and y coordinates.
    """
    x_start_adjusted = x_start - transition_dist1 * (1 - line_stop_frac)
    x_end_adjusted = x_end + transition_dist2 * line_stop_frac
    x_points = np.arange(x_start_adjusted, x_end_adjusted + resolution, resolution)

    x_middle = (x_start + x_end) / 2
    y_diff = y_end - y_start
    y_middle = y_start + (y_diff / 2)

    y_start_adjusted = y_middle - (y_diff / 2) * gradient / (x_middle - x_start_adjusted)
    y_end_adjusted = y_middle + (y_diff / 2) * gradient / (x_middle - x_end_adjusted)

    y_points = np.interp(x_points, [x_start_adjusted, x_end_adjusted], [y_start_adjusted, y_end_adjusted])

    return x_points, y_points

contours = gpd.read_file("alpine_contours_truncated_south.geojson")
end_points = gpd.read_file(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\howarth_alpine\endpoints\Alpine southern combined_truncated_dip_endpoints.shp")

min_transition_dist = 4000.  # meters
max_transition_dist = 10.e3  # meters
transition_dist_proportion = 0.25
spline_res = 20.  # meters
out_spline_res = 500.  # meters
spline_gradient = 1.2  # gradient of the spline
spline_k = 3  # spline order
spline_s = 0.0  # spline smoothing factor

end_points.rake[end_points.rake < 0.] += 360.  # Ensure rake values are non-negative
depths_to_match = np.intersect1d(contours['depth'], end_points['depth'])
trace = contours[contours['depth'] == 0.].geometry.values[0]
trace_end_points = end_points[end_points['depth'] == 0.]
trace_end_points["proj_distance"] = [trace.project(point) for point in trace_end_points.geometry]

transition_distances_dict = {}
trace_seg_lengths = np.array(trace_end_points.proj_distance)[1::2] - np.array(trace_end_points.proj_distance)[0:-1:2]
for i, seg_length in enumerate(trace_seg_lengths):
    default_transition_dist = transition_dist_proportion * seg_length
    if default_transition_dist < min_transition_dist:
        transition_dist = min_transition_dist
    elif default_transition_dist > max_transition_dist:
        transition_dist = max_transition_dist
    else:
        transition_dist = default_transition_dist

    transition_distances_dict[i] = transition_dist


step_start_dicts = {}
step_end_dicts = {}

output_list = []

for depth in depths_to_match:
    matched_contour = contours[contours['depth'] == depth].geometry.values[0]
    matched_end_points = end_points[end_points['depth'] == depth]
    matched_end_points["proj_distance"] = [matched_contour.project(point) for point in matched_end_points.geometry]
    step_start_dict = {}
    for i, dist in enumerate(matched_end_points.proj_distance):
        if i % 2 == 0:
            step_start_dict[dist] = int(i/2)

    step_end_dict = {}
    for i, dist in enumerate(matched_end_points.proj_distance):
        if i % 2 == 1:
            step_end_dict[dist] = int((i - 1)/2)
    step_start_dicts[depth] = step_start_dict
    step_end_dicts[depth] = step_end_dict

for depth in depths_to_match:
    matched_contour = contours[contours['depth'] == depth].geometry.values[0]
    matched_end_points = end_points[end_points['depth'] == depth]
    matched_end_points["proj_distance"] = [matched_contour.project(point) for point in matched_end_points.geometry]
    spline_x_list_sr = []
    spline_y_list_sr = []
    spline_x_list_rake = []
    spline_y_list_rake = []
    step_start_dict = step_start_dicts[depth]
    step_end_dict = step_end_dicts[depth]

    start_sr = matched_end_points.slip_rate.iloc[0]
    start_rake = matched_end_points.rake.iloc[0]
    start_transition_dist = transition_distances_dict[0]
    spline_x_sr, spline_y_sr = spline_points_start(matched_end_points.proj_distance.iloc[0], start_transition_dist, y_scaling=start_sr,
                                                    line_stop_frac=0.5, before_gradient=spline_gradient, resolution=spline_res)
    spline_x_list_sr.append(spline_x_sr)
    spline_y_list_sr.append(spline_y_sr)

    for i in range(len(matched_end_points) - 1):
        start_point = matched_end_points.geometry.iloc[i]
        end_point = matched_end_points.geometry.iloc[i + 1]
        start_distance = matched_end_points.proj_distance.iloc[i]
        end_distance = matched_end_points.proj_distance.iloc[i + 1]
        start_sr = matched_end_points.slip_rate.iloc[i]
        end_sr = matched_end_points.slip_rate.iloc[i + 1]
        start_rake = matched_end_points.rake.iloc[i]
        end_rake = matched_end_points.rake.iloc[i + 1]
        
        seg_length = end_distance - start_distance
        step_bool = any([start_sr != end_sr, start_rake != end_rake])
        if step_bool:
            start_transition_dist = transition_distances_dict[step_end_dict[start_distance]]
            end_transition_dist = transition_distances_dict[step_start_dict[end_distance]]
            spline_x_sr, spline_y_sr = spline_points_middle(start_distance, end_distance, start_transition_dist, end_transition_dist,
                                                    start_sr, end_sr, resolution=spline_res, gradient=spline_gradient)
            spline_x_rake, spline_y_rake = spline_points_middle(start_distance, end_distance, start_transition_dist, end_transition_dist,
                                                    start_rake, end_rake, resolution=spline_res, gradient=spline_gradient)
            spline_x_list_rake.append(spline_x_rake)
            spline_y_list_rake.append(spline_y_rake)
            spline_x_list_sr.append(spline_x_sr)
            spline_y_list_sr.append(spline_y_sr)
            if any(np.isnan(spline_y_sr)):
                # Handle NaN values in spline_y_sr
                print(f"NaN values found in slip rate spline for depth {depth} between distances {start_distance} and {end_distance}.")
        else:
            if start_distance not in step_start_dict.keys() or end_distance not in step_end_dict.keys():
                spline_x_rake = np.arange(end_distance, start_distance, spline_res)
                spline_y_rake = np.ones_like(spline_x_rake) * start_rake
                spline_x_sr = np.arange(end_distance, start_distance, spline_res)
                spline_y_sr = np.ones_like(spline_x_sr) * start_sr
            else:
                start_transition_dist = transition_distances_dict[step_start_dict[start_distance]]
                end_transition_dist = transition_distances_dict[step_end_dict[end_distance]]
                spline_x_rake = np.arange(start_distance + start_transition_dist, end_distance - end_transition_dist + spline_res, spline_res)
                spline_y_rake = np.ones_like(spline_x_rake) * start_rake
                spline_x_sr = np.arange(start_distance + start_transition_dist, end_distance - end_transition_dist + spline_res, spline_res)
                spline_y_sr = np.ones_like(spline_x_sr) * start_sr
            spline_x_list_rake.append(spline_x_rake)
            spline_y_list_rake.append(spline_y_rake)
            spline_x_list_sr.append(spline_x_sr)
            spline_y_list_sr.append(spline_y_sr)

    end_sr = matched_end_points.slip_rate.iloc[-1]
    end_rake = matched_end_points.rake.iloc[-1]
    end_transition_dist = transition_distances_dict[step_end_dict[matched_end_points.proj_distance.iloc[-1]]]
    spline_x_sr, spline_y_sr = spline_points_end(matched_end_points.proj_distance.iloc[-1], end_transition_dist, y_scaling=end_sr,
                                                    line_stop_frac=0.5, after_gradient=spline_gradient, resolution=spline_res)  
    spline_x_list_sr.append(spline_x_sr)
    spline_y_list_sr.append(spline_y_sr)

    spline_x_rake = np.concatenate(spline_x_list_rake)
    spline_y_rake = np.concatenate(spline_y_list_rake)
    spline_x_sr = np.concatenate(spline_x_list_sr)
    spline_y_sr = np.concatenate(spline_y_list_sr)

    spline_x_out = np.arange(matched_end_points.proj_distance.min(), matched_end_points.proj_distance.max() + out_spline_res, out_spline_res)
    spline_y_out_rake = make_interp_spline(spline_x_rake, spline_y_rake, k=spline_k)(spline_x_out)
    spline_y_out_sr = make_interp_spline(spline_x_sr, spline_y_sr, k=spline_k)(spline_x_out)
    spline_y_out_sr[spline_y_out_sr < 0] = 0.
    spline_y_out_rake[spline_y_out_rake > 180] -= 360.

    out_points = np.vstack([point.coords for point in matched_contour.interpolate(spline_x_out)])
    out_points = np.column_stack([out_points, spline_y_out_sr, spline_y_out_rake])
    output_list.append(out_points)



    # import matplotlib.pyplot as plt
    # plt.close('all')
    # plt.figure(figsize=(10, 6))
    # # plt.plot(spline_x_out, spline_y_out_rake, label='Rake', color='blue')
    # plt.plot(spline_x_out, spline_y_out_sr, label='Slip Rate', color='red')
    # plt.xlabel('Projected Distance')
    # plt.ylabel('Value')
    # plt.title('Interpolated Spline for Rake and Slip Rate')
    # # plt.legend()
    # plt.grid()
    # plt.show()

out_arr = np.vstack(output_list)
np.save("matched_slip_contours_truncated_south.npy", out_arr)
