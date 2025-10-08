import numpy as np
from scipy.spatial import KDTree
from rsqsim_api.fault.multifault import RsqSimMultiFault

# Read in the fault file
trimmed_faults = RsqSimMultiFault.read_fault_file_keith(fault_file="no_dip_change.flt")
# Get centres of all triangle patches
all_centres = np.vstack([fault.get_patch_centres() for fault in trimmed_faults.faults])

# Find nearby patches using KDTrees
centre_tree = KDTree(all_centres)
results = centre_tree.query_ball_point(all_centres, 7500.)
# First column of neighbours is the number of neighbours for each patch
num_neighbours = [len(item) for item in results]
# Write out neighbour file
with open("neighbour_file_no_dip_change.txt", "w") as f:
    for i, row in enumerate(results):
        for item in row[:-1]:
            if item != i:
                f.write(f"{item:d} ")
        f.write(f"{row[-1]:d}\n")