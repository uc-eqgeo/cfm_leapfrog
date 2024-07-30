import numpy as np

southern_xyz = np.array([1530472.4, 5282209.2, 0.])
transition_xyz = np.array([2.26287e6, 5.90386e6, 0.])
northern_xyz = np.array([2688695.6, 7299458.5, 0.])

strike = 26.
strike_vector = np.array([np.sin(np.deg2rad(strike)), np.cos(np.deg2rad(strike)), 0.])


box_width = 5.e5
buffer = 1.e5

def box_params(start_centre, end_centre, strike_vector, box_width, buffer, box_height=2.e5):

    start_end_dist_along_strike = np.dot(end_centre - start_centre, strike_vector)
    centre_dist = start_end_dist_along_strike / 2.
    print(centre_dist)
    if start_end_dist_along_strike < 0:
        centre_dist -= buffer/2.
    else:
        centre_dist += buffer/2.

    centre = start_centre + centre_dist * strike_vector
    length = np.abs(start_end_dist_along_strike) + buffer
    print("Centre:", centre)
    print("Length:", box_width, length, box_height)
    print("x-axis:", strike_vector[1], -strike_vector[0], 0.)
    print("y-axis:", strike_vector[0], strike_vector[1], 0.)
    print("z-axis:", 0., 0., 1.)

print("Southern cutting box")
box_params(southern_xyz, southern_xyz - np.array([0., 1., 0.]) * 5.e4, np.array([0., 1., 0.]), box_width, buffer)

print("Centre cutting box")
box_params(transition_xyz, transition_xyz + 150e3 * strike_vector, strike_vector, box_width, 0.)

print("Southern mid box 1")
box_params(transition_xyz, southern_xyz, strike_vector, 1.e6, buffer)

print("Northern mid box 1")
box_params(transition_xyz + 50e3 * strike_vector, northern_xyz, strike_vector, 1.e6, buffer)

print("Southern mid box 2")
box_params(transition_xyz + 50e3 * strike_vector, southern_xyz, strike_vector, 1.e6, buffer)

print("Northern mid box 2")
box_params(transition_xyz + 100e3 * strike_vector, northern_xyz, strike_vector, 1.e6, buffer)

print("Southern mid box 3")
box_params(transition_xyz + 100e3 * strike_vector, southern_xyz, strike_vector, 1.e6, buffer)

print("Northern mid box 3")
box_params(transition_xyz + 150e3 * strike_vector, northern_xyz, strike_vector, 1.e6, buffer)

