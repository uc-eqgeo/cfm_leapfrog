import igl
import numpy as np
import scipy
import glob
import os.path

min_area = 1.e5
stls = ["Puysegur_second_trimmed_lt40_remeshed_2500.stl"]
epsilons = np.array([10**x for x in range(-6, -2)])

for stl in stls:
    v, f = igl.read_triangle_mesh(stl)
    min_tri_area = (igl.doublearea(v, f) / 2.0).min()

    if min_tri_area < min_area:
        for eps in epsilons:
            out = igl.collapse_small_triangles(v, f, eps=eps)
            min_tri_area_i = (igl.doublearea(v, out) / 2.0).min()
            if min_tri_area_i > min_area:
                break
        igl.write_triangle_mesh(os.path.basename(stl).replace(".stl", "_collapsed.stl"), v, out)
        print(stl)
        print(f"Collapsed {len(f) - len(out)} triangles")
    else:
        igl.write_triangle_mesh(os.path.basename(stl).replace(".stl", "_collapsed.stl"), v, f)


