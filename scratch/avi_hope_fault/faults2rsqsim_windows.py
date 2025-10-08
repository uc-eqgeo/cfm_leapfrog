# Import relevant modules
import matplotlib.pyplot as plt
import shapely.plotting
from rsqsim_api.fault.multifault import RsqSimMultiFault, RsqSimSegment
from rsqsim_api.io.rsqsim_constants import seconds_per_year
import os, glob
import pandas as pd
import geopandas as gpd
import numpy as np
import difflib
import csv
import shapely as sp
from shapely.geometry import LineString,Point
from shapely.ops import unary_union
from shapely.plotting import plot_points
from scipy.interpolate import griddata,UnivariateSpline
from scipy.spatial import KDTree
import meshio

def moving_average(a,n=10):
    """calculate moving average over specified number of points"""
    cumsum=np.cumsum(a, dtype=float)
    cumsum[n:] =cumsum[n:] - cumsum[:-n]
    return cumsum[n-1:]/n

def make_spline(val1,dist1,val2: float =-1.,dist2: float = 0.,grad_factor: float = 0.8):
    total_dist = dist1 + dist2
    #spline should end up at average value if 2 are given or default is to taper to 0
    if val2>=0.:
        val2taper = (val1+val2)/2.
    else:
        val2taper = 0.
    x_points_near_centre = np.linspace(0, dist1 / 20., num=10)
    grad = (val1 - val2taper) * grad_factor / dist1
    centre_vals = x_points_near_centre * grad + val2taper
    starter_x_points = np.arange(dist1, 1.2 * dist1, 20)
    added_y_points = np.ones(len(starter_x_points)) * val1
    x_points = np.append(x_points_near_centre, starter_x_points)
    y_vals = np.append(centre_vals, added_y_points)
    spline = UnivariateSpline(x_points, y_vals, k=3, s=0)
    return spline

# Tell python where field paths etc are relative to
script_dir = os.path.abspath('')
#NB vtks will be written to this directory
mesh_dir =r"remeshed_faults"
#directory to write rsqsim fault file to
out_dir=r"rsqsim_inputs_new_rake"
# read in cfm to get slip rates
cfm = gpd.read_file(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\avi_hope_fault\hope_and_alpine.geojson")
#proportion of downdip width to start tapering slip from for a) base and b) fault edges
base_thresh_factor=1.
side_thresh_factor=0.25
min_side_dist = 4000. #maximum distance from fault edges to taper slip over
max_side_dist = 10000. #maximum distance from fault edges to taper slip over
#increment along fault to check patch locations
trace_spacing = 2000.
#read in fault combination file to separate back out

with open(r"C:\Users\arh128\PycharmProjects\cfm_leapfrog\scratch\avi_hope_fault\hope_alpine_connected_edited.csv") as fin:
    fault_combs=csv.reader(fin)
    fault_comb_dict={}
    for line in fault_combs:
        fault_comb=line[0]
        faults=line[1:]
        fault_comb_dict[fault_comb]=faults

all_faults = []
#problem_faults=['Chancet','Akatarawa','CampbellBank','Papatea']
small_patch_faults =[]
for i, file in enumerate(glob.iglob(f'{mesh_dir}/*.stl')):
#for i, file in enumerate(glob.iglob(f'{mesh_dir}/BigBen_remeshed.stl')):
    fltName = os.path.basename(file).split('.')[0]
    fltName = fltName.replace(" ","_")
    segNum = i + 1
    print(fltName)
    # #if segNum > 100:
    if "combined" in fltName:

        #combined faults need to be re-separated to match up cfm slip rates
        # create fault segment first
        basename = difflib.get_close_matches(fltName, fault_comb_dict.keys(), n=1)[0]
        #try:
        fltseg = RsqSimSegment.from_stl(stl_file=file, segment_number=segNum, fault_name=fltName,
                                        strike_slip=0., dip_slip=0.)
        edge_patches_nos = fltseg.find_edge_patch_numbers(top=False)
        wholeTrace = LineString(fltseg.find_top_vertices(depth_tolerance=10.))
        wholeFault = {'Name': [fltName], 'geometry': [wholeTrace]}
        wholegdf = gpd.GeoDataFrame.from_dict(wholeFault, geometry=wholeFault['geometry'], crs="EPSG:2193")

        fltNames = fault_comb_dict[basename]
        fltDict = {"Name": [], "geometry": [], "slip rate": [], "rake": []}
        for flt in fltNames:
            if len(flt) > 0:
                #find closest cfm name to fault section
                newFltName = difflib.get_close_matches(flt, cfm['Name'])[0]
                fltDict["Name"].append(flt)
                fltDict["geometry"].append(cfm['geometry'].loc[cfm['Name'] == newFltName].values[0])
                fltDict["slip rate"].append(cfm['SR_pref'].loc[cfm['Name'] == newFltName].values[0])
                fltDict['rake'].append(cfm['Rake_pref'].loc[cfm['Name'] == newFltName].values[0])

        cfmSegs = gpd.GeoDataFrame.from_dict(fltDict, geometry=fltDict["geometry"], crs="EPSG:2193")
        #find spatial join between cfm and trace based on nearest neighbours
        intersecting_segs = cfmSegs.sjoin_nearest(wholegdf)
        # make set up of evenly spaced points along each of these traces

        cfmSegs["trace_dists"]=[np.arange(0.,sp.length(line),trace_spacing) for line in cfmSegs['geometry']]
        cfmSegs["trace_points"]=[cfmSegs["geometry"].iloc[i].interpolate(dist) for i,dist in enumerate(cfmSegs["trace_dists"]) ]



        # calculate distance from end of each segment trace where tapering should start
        taper_dists=[]
        for line in cfmSegs['geometry']:
            seg_length=sp.length(line)
            seg_taper_potential = seg_length * side_thresh_factor
            if seg_taper_potential < max_side_dist:
                if seg_taper_potential > min_side_dist:
                    side_taper_dist = seg_length * side_thresh_factor
                else:
                    side_taper_dist = min_side_dist
            else:
                side_taper_dist = max_side_dist
            taper_dists.append(side_taper_dist)
        cfmSegs["taper_dists"]=taper_dists
        #find adjacent segments + associated average slip
        cfmSegs_buffer = cfmSegs.copy()
        #add buffer to allow for slight discrepancies in trace work
        cfmSegs_buffer.geometry = cfmSegs.geometry.buffer(200)
        neighbourCFMSegs=cfmSegs_buffer.sjoin(cfmSegs_buffer,predicate="intersects",how='left')
        #neighbourCFMSegs.reindex(np.arange(0,len(neighbourCFMSegs)),axis=0)
        #average slip at each intersection
        neighbourCFMSegs["average_slip"]=neighbourCFMSegs[['slip rate_left','slip rate_right']].mean(axis=1)
        #distance over which transition occurs
        neighbourCFMSegs["edge_distance"] = neighbourCFMSegs[['taper_dists_left','taper_dists_right']].sum(axis=1)
        # set of points for spline to fit through
        splines=[]
        for idx in np.arange(0,len(neighbourCFMSegs),1):
            slip1=neighbourCFMSegs["slip rate_left"].iloc[idx]
            # for sideways splines
            taper_dist=neighbourCFMSegs["taper_dists_left"].iloc[idx]
            average_slip=neighbourCFMSegs["average_slip"].iloc[idx]
            total_dist=neighbourCFMSegs["edge_distance"].iloc[idx]
            x_points_near_centre = np.linspace(0,total_dist/20.,num=10)
            grad = (slip1-average_slip)*1.1 /taper_dist
            slip_vals_centre= x_points_near_centre*grad + average_slip
            starter_x_points = np.arange(taper_dist, 1.2*taper_dist,20)
            added_y_points = np.ones(len(starter_x_points))*slip1
            x_points = np.append(x_points_near_centre,starter_x_points)
            slip_vals = np.append(slip_vals_centre,added_y_points)
            spline=UnivariateSpline(x_points,slip_vals,k=3,s=0)
            splines.append(spline)


        neighbourCFMSegs["splines"]=splines

        # make a base spline from 1 to 0 over the threshold distance to multiply slip rates by


        # now make splines for either end of the whole fault (tapering to 0 over max side dist)
        # end segments only appear once in the Names_left column of neighbouringCFMSegs

        # end_spline_dict = {}
        end_names=neighbourCFMSegs["Name_left"].value_counts().keys()[-2:]
        # for end_name in end_names:
        #     cfm_idx = np.where(cfmSegs["Name"]==end_name)
        #     taper_dist = cfmSegs["taper_dists"].iloc[cfm_idx].values[0]
        #     slip_val = cfmSegs["slip rate"] .iloc[cfm_idx].values[0]
        #     x_points_near_centre = np.linspace(0, taper_dist / 20., num=10)
        #     grad = (slip_val) * 1.1 / taper_dist
        #     slip_vals_centre = x_points_near_centre * grad
        #     starter_x_points = np.arange(taper_dist, 1.2 * taper_dist, 20)
        #     added_y_points = np.ones(len(starter_x_points)) * slip_val
        #     x_points = np.append(x_points_near_centre, starter_x_points)
        #     slip_vals = np.append(slip_vals_centre, added_y_points)
        #     spline = UnivariateSpline(x_points, slip_vals, k=3, s=0)
        #     end_spline_dict[end_name]=spline


        #now find values on actual patches


        patches = fltseg.patch_dic.values()
        patch_ids = list(fltseg.patch_dic.keys())
        #find patch centres
        centres=[Point(patch.centre) for patch in patches]
        # extra stuff for JKN
        # if fltName == "JordanKekerenguNeedlescombined":
        #     print("JKN special case")
        #     edge_patches_nos = fltseg.find_edge_patch_numbers(top=False)
        #
        #     no_geom_centres = [patch.centre for patch in patches]
        #     edge_centres = [patch.centre for i, patch in enumerate(patches) if patch_ids[i] in edge_patches_nos]
        #
        #     # create KDTree for efficient nearest neighbour look up
        #     edge_tree = KDTree(edge_centres)
        #     # distances to edge
        #     dists2edge = edge_tree.query(no_geom_centres)[0]


        #top edge of
        top_edge = fltseg.trace
        #create points at 1km spacing along trace
        trace_dists=np.arange(0.,sp.length(top_edge),trace_spacing)
        trace_points=[top_edge.interpolate(dist) for dist in trace_dists]

        # make array of [dist along whole trace, distance along seg trace, depth, proportion of depth] for each patch
        patch_loc_in_trace_coords=np.zeros([len(patches),4])
        #patch_neighbours=np.zeros([len(patches),2])
        #first look at whole trace
        for i,patch in enumerate(patches):

            #find distance along whole trace
            distances2trace=centres[i].distance(trace_points)
            along_trace_dist=trace_dists[np.where(distances2trace==np.min(distances2trace))]

            patch_loc_in_trace_coords[i, 0] = along_trace_dist
            patch_loc_in_trace_coords[i, 2] = -1*patch.centre[2]


        # # find depth at each trace dist
        # dist_depth = np.zeros([len(trace_dists), 2])
        # for i, dist in enumerate(trace_dists):
        #     dist_depth[i, 0] = dist
        #     if any(patch_loc_in_trace_coords[:, 0] == dist):
        #         dist_depth[i, 1] = np.max(
        #             patch_loc_in_trace_coords[np.where(patch_loc_in_trace_coords[:, 0] == dist), 2])
        # smooth to get rid of zeros and make a "base trace"
        # non_zero_depths = dist_depth[np.where(dist_depth[:, 1] != 0)]
        # # take moving average
        # if len(non_zero_depths) > 100.:
        #     n = 10
        #     non_zero_depths_smooth = moving_average(non_zero_depths[:, 1], n=n)
        #     depth_alongstrike_dists = non_zero_depths[n - 1:, 0]
        #     depth_trace_spline = UnivariateSpline(depth_alongstrike_dists, non_zero_depths_smooth, ext=3, s=0)
        # elif len(non_zero_depths) > 10:
        #     n = 3
        #     non_zero_depths_smooth = moving_average(non_zero_depths[:, 1], n=n)
        #     depth_alongstrike_dists = non_zero_depths[n - 1:, 0]
        #     depth_trace_spline = UnivariateSpline(depth_alongstrike_dists, non_zero_depths_smooth, ext=3, s=0)
        # else:
        #     depth_trace_spline = UnivariateSpline(non_zero_depths[:, 0], non_zero_depths[:, 1])

        centre_slip=[]
        for i,patch in enumerate(patches):
            # matching up to CFM
            cfm_dists = [centres[i].distance(line) for line in cfmSegs['geometry']]
            min_cfm_dist = np.min(cfm_dists)
            idx = np.where(cfm_dists == min_cfm_dist)[0]
            if len(idx) > 1:
                idx = [idx[0]]

            # find name of closest segment and next closest segment
            close_seg = cfmSegs['Name'].iloc[idx].values[0]

            total_slip = cfmSegs['slip rate'].iloc[idx].values[0]
            rake = cfmSegs['rake'].iloc[idx].values[0]
            seg_length = sp.length(cfmSegs['geometry'].iloc[idx].values[0])

            seg_trace_points = cfmSegs['trace_points'].iloc[idx].values[0]
            seg_trace_dists = cfmSegs['trace_dists'].iloc[idx].values[0]

            #find distance to each segment then find the minimum
            distances2segtrace = centres[i].distance(seg_trace_points)
            along_seg_dist = seg_trace_dists[np.where(distances2segtrace == np.min(distances2segtrace))]
            DistFromSegEnd = np.min([along_seg_dist, seg_length - along_seg_dist])
            #find distance to trace end
            DistFromTraceEnd = np.min([patch_loc_in_trace_coords[i, 0],(sp.length(wholeTrace) - patch_loc_in_trace_coords[i, 0])])




            patch_loc_in_trace_coords[i, 1] = along_seg_dist

            # #and find distance to base of fault at that location
            # local_max_depth = depth_trace_spline(patch_loc_in_trace_coords[i,0])
            # base_thresh_dist = local_max_depth*base_thresh_factor
            # DistFromBase = local_max_depth - (-1.* patch.centre[2])
            # if DistFromBase < 0:
            #     DistFromBase = 0.

            #patch_neighbours[i, 0] = DistFromSegEnd
            #patch_neighbours[i, 1] = DistFromBase
            # now assign slip
            #is the patch on an end segment and close to the end of the fault?

            # if close_seg == "Jordan" and DistFromTraceEnd <= seg_length/2.:
            #     #print("Special case for Jordan")
            #     boundary_spline = end_spline_dict[close_seg]
            #     edge_dist = dists2edge[i]
            #     slip_rate = boundary_spline(edge_dist)
            #elif
            if close_seg in end_names and DistFromTraceEnd <side_taper_dist:
                # boundary_spline = end_spline_dict[close_seg]
                # slip_rate = boundary_spline(DistFromTraceEnd)
                slip_rate = total_slip

            elif DistFromSegEnd < side_taper_dist and seg_length > 2*max_side_dist:
                # is the patch close to the end of a segment (which isn't too short to taper across)?
                #taper to side
                #find next closest segment
                if len(cfmSegs) == 2:
                    next_close_seg_idx = [i for i in cfmSegs.index if cfmSegs['Name'].iloc[i] != close_seg]
                else :
                    next_close_seg_idx = np.where(cfm_dists == np.partition(cfm_dists, 2)[1])[0]
                    if len(next_close_seg_idx) > 1:
                        # if there are two traces equidistant from the patch use the one which hasnt already been assigned as closest (not very robust!)
                        next_close_seg_idx = [next_close_seg_idx[1]]
                next_close_seg = cfmSegs['Name'].iloc[next_close_seg_idx].values[0]
                spline_idx=np.where(np.logical_and(neighbourCFMSegs['Name_left']==close_seg, neighbourCFMSegs['Name_right'] == next_close_seg))
                boundary_spline=neighbourCFMSegs["splines"].iloc[spline_idx].values[0]

                slip_rate = boundary_spline(DistFromSegEnd)
            else:
                slip_rate = total_slip

            # if DistFromBase < base_thresh_dist: # and not (close_seg == "Jordan" and DistFromTraceEnd <= seg_length/2.):
            #     base_spline = make_spline(1,base_thresh_dist)
            #     slip_rate=slip_rate*base_spline(DistFromBase)
            if slip_rate < 1e-3:
                slip_rate = 0.
            # # check no patches at edge have full slip rate
            # if patch_ids[i] in edge_patches_nos:
            #     if slip_rate > total_slip / 2.:
            #         slip_rate = slip_rate / 2.


            patch.strike_slip = slip_rate * np.cos(np.radians(rake))
            patch.dip_slip = slip_rate * np.sin(np.radians(rake))

            centre_slip.append(np.hstack([patch.centre, slip_rate]))

        # mesh = meshio.Mesh(points=fltseg.vertices, cells=[("triangle", fltseg.triangles)])
        # mesh.cell_data["distFromSegEnd"] = patch_neighbours[:,0]
        # mesh.cell_data["distFrombase"] = patch_neighbours[:, 1]
        # mesh.write(os.path.join(mesh_dir, f'{fltName}_neighbours.vtk'),file_format='vtk')
        fltseg.to_vtk(os.path.join(mesh_dir, f'{fltName}.vtk'), write_slip=True)
        fltseg.to_rsqsim_fault_file(os.path.join(mesh_dir, f'{fltName}.flt'))
        #write slip in rbf interpolant friendly format
        flt_slips = np.array(centre_slip, dtype='float')
        np.savetxt(os.path.join(mesh_dir, f'{fltName}_patch_slips.csv'), flt_slips, fmt=['%.6e', '%.6e', '%.6e', '%.1f'],
                   delimiter=',', newline='\n', header='E N X slip_rate')
        fltarray = fltseg.to_rsqsim_fault_array()
        all_faults.append(fltarray)
        #
        #
        patch_areas = [patch.area for patch in patches]
        if any([area <1.e4 for area in patch_areas]):
            print(f'Fault {basename} has small patches.')
            small_patch_faults.append(fltName)

    # elif fltName in problem_faults:
    #     print(f"{fltName} in problem faults")
    #     #use distance to nearest edge patch as distance for faults which don't reach the surface everywhere
    #     newFltName = difflib.get_close_matches(fltName, cfm['Name'])[0]
    #     fltseg = RsqSimSegment.from_stl(stl_file=file, segment_number=segNum, fault_name=fltName,
    #                                     strike_slip=0., dip_slip=0.)
    #     total_slip = cfm['SR_pref'].loc[cfm['Name'] == newFltName].values[0]
    #     #print(total_slip)
    #     rake = float(cfm['Rake_pref'].loc[cfm['Name'] == newFltName].values)
    #
    #     patches = fltseg.patch_dic.values()
    #     patch_ids = list(fltseg.patch_dic.keys())
    #     edge_patches_nos = fltseg.find_edge_patch_numbers(top=False)
    #
    #     centres = [patch.centre for patch in patches]
    #     edge_centres = [patch.centre for i,patch in enumerate(patches) if patch_ids[i] in edge_patches_nos]
    #
    #     #create KDTree for efficient nearest neighbour look up
    #     edge_tree = KDTree(edge_centres)
    #     #distances to edge
    #     dists2edge = edge_tree.query(centres)[0]
    #
    #
    #     #find side taper distance
    #     top_edge = fltseg.trace
    #     trace_length = sp.length(top_edge)
    #
    #     #make side spline
    #     if trace_length > max_side_dist:
    #
    #         seg_taper_potential = trace_length * side_thresh_factor
    #         if seg_taper_potential < max_side_dist:
    #             if seg_taper_potential > min_side_dist:
    #                 side_taper_dist = seg_taper_potential
    #             else:
    #                 side_taper_dist = min_side_dist
    #         else:
    #             side_taper_dist = max_side_dist
    #
    #         side_spline=make_spline(total_slip,side_taper_dist)
    #
    #         #calculate slip rate based on distance to edge patch
    #         slip_rates = [side_spline(dist2edge) for dist2edge in dists2edge]
    #         centre_slip=[]
    #         for i,patch in enumerate(patches):
    #             slip_rate = slip_rates[i]
    #             if slip_rate < 1e-3:
    #                 slip_rate = 0.
    #             patch.rake = rake
    #             patch.strike_slip = slip_rate * np.cos(np.radians(patch.rake))
    #             patch.dip_slip = slip_rate * np.sin(np.radians(patch.rake))
    #
    #             centre_slip.append(np.hstack([patch.centre, slip_rate]))
    #
    #     else:
    #         #if have a short fault, just set outer patches (except top) to half the slip rate
    #         centre_slip=[]
    #         for patch_id in fltseg.patch_dic.keys():
    #             patch = fltseg.patch_dic[patch_id]
    #             if patch_id in edge_patches_nos:
    #                 slip_rate = total_slip/2.
    #             else:
    #                 slip_rate = total_slip
    #             if slip_rate < 1e-3:
    #                 slip_rate = 0.
    #             patch.rake = rake
    #             patch.strike_slip = slip_rate * np.cos(np.radians(patch.rake))
    #             patch.dip_slip = slip_rate * np.sin(np.radians(patch.rake))
    #
    #             centre_slip.append(np.hstack([patch.centre, slip_rate]))
    #
    #     patch_areas = [patch.area for patch in patches]
    #     if any([area < 1.e4 for area in patch_areas]):
    #         print(f'Fault {fltName} has small patches')
    #         small_patch_faults.append(fltName)
    #
    #     fltseg.to_vtk(os.path.join(mesh_dir, f'{fltName}.vtk'), write_slip=True)
    #     fltseg.to_rsqsim_fault_file(os.path.join(mesh_dir, f'{fltName}.flt'))
    #     fltarray = fltseg.to_rsqsim_fault_array()
    #     all_faults.append(fltarray)
    #     # write slip in rbf interpolant friendly format
    #     flt_slips = np.array(centre_slip, dtype='float')
    #     np.savetxt(os.path.join(mesh_dir, f'{fltName}_patch_slips.csv'), flt_slips,
    #                fmt=['%.6e', '%.6e', '%.6e', '%.1f'],
    #                delimiter=',', newline='\n', header='E N X slip_rate')

    else:




        newFltName = difflib.get_close_matches(fltName, cfm['Name'])[0]
        fltseg = RsqSimSegment.from_stl(stl_file=file, segment_number=segNum, fault_name=fltName,
                                        strike_slip=0., dip_slip=0.)
        total_slip = cfm['SR_pref'].loc[cfm['Name'] == newFltName].values[0]

        rake = float(cfm['Rake_pref'].loc[cfm['Name'] == newFltName].values)

        patches = fltseg.patch_dic.values()
        patch_ids = list(fltseg.patch_dic.keys())
        edge_patches_nos = fltseg.find_edge_patch_numbers(top=False) #edge_patches_nos = fltseg.find_edge_patch_numbers(top=False)

        centres = [Point(patch.centre) for patch in patches]

        # #find side taper distance
        # top_edge = fltseg.trace
        # trace_length = sp.length(top_edge)
        #
        # if trace_length > max_side_dist:
        #
        #     seg_taper_potential = trace_length * side_thresh_factor
        #     if seg_taper_potential < max_side_dist:
        #         if seg_taper_potential > min_side_dist:
        #             side_taper_dist = seg_taper_potential
        #         else:
        #             side_taper_dist = min_side_dist
        #     else:
        #         side_taper_dist = max_side_dist
        #
        #     side_spline=make_spline(total_slip,side_taper_dist)
        #
        #     # create points at regular spacing along trace
        #     trace_dists = np.arange(0., sp.length(top_edge), trace_spacing)
        #     trace_points = [top_edge.interpolate(dist) for dist in trace_dists]
        #
        #     patch_loc_in_trace_coords = np.zeros([len(patches), 3])
        #     # patch_neighbours = np.zeros([len(patches), 2])
        #     slip_rates = []
        #     # first look at patch depths along trace
        #     for i, patch in enumerate(patches):
        #         # find distance along whole trace
        #         distances2trace = centres[i].distance(trace_points)
        #         along_trace_dist = trace_dists[np.where(distances2trace == np.min(distances2trace))]
        #         #print(along_trace_dist)
        #         DistFromSegEnd = np.min([along_trace_dist,trace_length - along_trace_dist])
        #         #print(DistFromSegEnd)
        #         # taper at sides
        #         if DistFromSegEnd < side_taper_dist:
        #             slip_rate = side_spline(DistFromSegEnd)
        #         else:
        #             slip_rate = total_slip
        #         if slip_rate < 1e-3:
        #             slip_rate = 0.
        #
        #         slip_rates.append(slip_rate)
        #         patch_loc_in_trace_coords[i, 0] = along_trace_dist
        #         patch_loc_in_trace_coords[i, 1] = -1 * patch.centre[2]

        # # find depth at each trace dist
        # dist_depth = np.zeros([len(trace_dists), 2])
        # for i, dist in enumerate(trace_dists):
        #     dist_depth[i, 0] = dist
        #     if any(patch_loc_in_trace_coords[:, 0] == dist):
        #         dist_depth[i, 1] = np.max(
        #             patch_loc_in_trace_coords[np.where(patch_loc_in_trace_coords[:, 0] == dist), 1])
        # # smooth to get rid of zeros and make a "base trace"
        # non_zero_depths = dist_depth[np.where(dist_depth[:, 1] != 0)]
        # # set uniform base depth if variation is small
        # if np.ptp(non_zero_depths) < 2500 or len(non_zero_depths < 10):
        #     base_depth = np.mean(non_zero_depths[:,1])
        # else:
        #     base_depth= 0.
        #     if len(non_zero_depths) > 50.:
        #         n = 10
        #         non_zero_depths_smooth = moving_average(non_zero_depths[:, 1], n=n)
        #         depth_alongstrike_dists = non_zero_depths[n - 1:, 0]
        #         depth_trace_spline = UnivariateSpline(depth_alongstrike_dists, non_zero_depths_smooth, ext=3, s=0)
        #     else:
        #         n = 3
        #         non_zero_depths_smooth = moving_average(non_zero_depths[:, 1], n=n)
        #         depth_alongstrike_dists = non_zero_depths[n - 1:, 0]
        #         depth_trace_spline = UnivariateSpline(depth_alongstrike_dists, non_zero_depths_smooth, ext=3, s=0)

        centre_slip=[]
        for i, patch in enumerate(patches):
            #slip_rate = slip_rates[i]
            slip_rate =total_slip
            # find distance to base of fault at that location
            # if base_depth >0.:
            #     local_max_depth = base_depth
            # else:
            #     local_max_depth = depth_trace_spline(patch_loc_in_trace_coords[i, 0])
            #
            # base_thresh_dist = local_max_depth * base_thresh_factor
            # DistFromBase = local_max_depth - (-1. * patch.centre[2])
            #
            # if DistFromBase < 0:
            #     DistFromBase = 0.
            # patch_neighbours[i, 1] = 0
            #taper at base
            # if DistFromBase < base_thresh_dist:
            #     base_spline = make_spline(1, base_thresh_dist)
            #     slip_rate = slip_rate * base_spline(DistFromBase)
                # patch_neighbours[i, 1] = base_spline(DistFromBase)
            if slip_rate < 1e-3:
                slip_rate = 0.

            #check no patches at edge have full slip rate
            # if patch_ids[i] in edge_patches_nos:
            #     # patch_neighbours[i,0] =1
            #     if slip_rate > total_slip/2.:
            #         slip_rate = slip_rate /2.

            patch.rake = rake
            patch.strike_slip = slip_rate * np.cos(np.radians(rake))
            patch.dip_slip = slip_rate * np.sin(np.radians(rake))

            centre_slip.append(np.hstack([patch.centre, slip_rate]))

                # mesh = meshio.Mesh(points=fltseg.vertices, cells=[("triangle", fltseg.triangles)])
                # mesh.cell_data["edge_patch"] = patch_neighbours[:,0]
                # mesh.cell_data["basedist"] = patch_neighbours[:, 1]
                # mesh.write(os.path.join(mesh_dir, f'{fltName}_neighbours.vtk'),file_format='vtk')
        # else:
        #     #if have a short fault, just set outer patches (except top) to half the slip rate
        #     centre_slip = []
        #     for patch_id in fltseg.patch_dic.keys():
        #         patch = fltseg.patch_dic[patch_id]
        #         # if patch_id in edge_patches_nos:
        #         #     slip_rate = total_slip/2.
        #         # else:
        #         #     slip_rate = total_slip
        #         slip_rate = total_slip
        #         patch.rake = rake
        #         patch.strike_slip = slip_rate * np.cos(np.radians(patch.rake))
        #         patch.dip_slip = slip_rate * np.sin(np.radians(patch.rake))
        #
        #         centre_slip.append(np.hstack([patch.centre, slip_rate]))

        fltseg.to_vtk(os.path.join(mesh_dir, f'{fltName}.vtk'), write_slip=True)
        fltseg.to_rsqsim_fault_file(os.path.join(mesh_dir, f'{fltName}.flt'))
        fltarray = fltseg.to_rsqsim_fault_array()
        all_faults.append(fltarray)
        # write slip in rbf interpolant friendly format
        flt_slips = np.array(centre_slip, dtype='float')
        np.savetxt(os.path.join(mesh_dir, f'{fltName}_patch_slips.csv'), flt_slips,
                   fmt=['%.6e', '%.6e', '%.6e', '%.1f'],
                   delimiter=',', newline='\n', header='E N X slip_rate')

        patch_areas = [patch.area for patch in patches]
        if any([area <1.e4 for area in patch_areas]):
            indices = np.where([area < 1.e4 for area in patch_areas])[0]
            print(f'Fault {fltName} has small patches.')
            small_patch_faults.append((fltName, indices))
            for i in indices:
                print(f'Patch {i} has area {patch_areas[i]}')



print(small_patch_faults)
fault_segs = pd.concat(all_faults)
fault_segs.to_csv(os.path.join(out_dir,"crustal_faults.flt"), index=False, header=False, sep='\t', encoding='ascii')

#if want to check outputs
# fault_model=RsqSimMultiFault.read_fault_file_keith(os.path.join(script_dir,'flt_files_v3/all_faults.flt'))
# fault_model.slip_rate_to_vtk(os.path.join(script_dir,"all_slip_rate.vtk"))