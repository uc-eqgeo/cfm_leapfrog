import numpy as np

west_slip_rate = 15.2 # mm/yr
west_strike = 82 # degrees
south_strike = 98 # degrees

south_dip = 70 # degrees

strike_diff = np.abs(west_strike - south_strike)
south_ss = west_slip_rate * np.cos(np.radians(strike_diff))
south_sp = west_slip_rate * np.sin(np.radians(strike_diff))
south_ds = south_sp * np.tan(np.radians(south_dip))

south_total = np.linalg.norm([south_ss, south_ds])