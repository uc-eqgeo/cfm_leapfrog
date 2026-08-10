"""Whole New Zealand fault-mesh pipeline.

Builds triangulated fault-surface meshes for the whole-of-New-Zealand
Community Fault Model, from raw fault traces to clean, uniformly-sized
triangle meshes. Three sequential stages:

1. Build surfaces - contour + triangulate each curated fault
                     -> test_objs/, test_vtks/, test_contours/
2. Cut surfaces    - trim against higher-priority faults and the base
                     depth surface -> test_final_meshes/*_cut.obj
2b. Check cuts     - flag cut meshes still reaching the surface well beyond
                     their original trace, and cuts that left a jagged
                     subsurface outline -> surface_trace_check.csv/.geojson,
                     jagged_cut_check.csv/.geojson
3. Remesh          - optimise every cut surface into uniform,
                     near-equilateral triangles with MMG
                     -> test_remeshed_vtks/, merged_mesh_remeshed.vtk

Run from this file's directory (scratch/whole_nz_no_leapfrog/) inside the
leapfrog-fault-models conda environment; paths below are relative to that
folder. Each stage reads its inputs back from disk, so a later stage can be
re-run on its own (e.g. re-cut without re-building the raw surfaces).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import meshio
import pyvista as pv

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.mesh import FaultMesh
from fault_mesh.io.array_operations import read_raster


def write_geojson_safe(gdf, path, **kwargs):
    """Write a GeoDataFrame to disk, skipping with a warning if the file is
    locked (e.g. still open as a QGIS layer) instead of crashing the run.
    These GeoJSONs are inspection-only outputs, so a skipped write just
    means the file keeps its previous contents until QGIS releases it."""
    try:
        gdf.to_file(str(path), **kwargs)
    except PermissionError as e:
        print(f"  Warning: could not write {path} (locked, e.g. open in QGIS?): {e}")

# --- Input files (relative to this script's folder) ---
FAULT_SHP         = Path("NZ_CFM_v1_0_rs1km_modified_GDM.gpkg")  # raw fault traces (whole-NZ CFM)
FAULT_SYSTEMS     = Path("NZ_CFM_v1_0_rs1km_modified_GDM_connected_edited.csv")  # which segments form one fault
CUTTING_HIERARCHY = Path("NZ_CFM_v1_0_rs1km_modified_GDM_hierarchy_edited.csv")  # which fault cuts which
DEPTH_RASTER      = Path("with_hannu_mods.tif")          # base/Moho-style depth surface (GeoTIFF, z in metres)
ADDITIONAL_CUTS   = Path("additional_cuts.csv")          # extra cut relationships not implied by the hierarchy
EXCLUDED_CUTS     = Path("excluded_cuts.csv")            # pairs that must NEVER be cut (optional; applied only if present)
TRACE_EXT_EDITED  = Path("edited_trace_extensions.csv")  # curated along-strike trace extensions (optional; applied only if present)
TRACE_EXT_REMOVED = Path("removed_trace_extensions.csv") # extensions ruled out in review (optional; subtracted from the above)

# --- Coordinate system & network-building parameters ---
EPSG             = 2193      # NZTM2000; set to None to skip reprojection
DIST_TOLERANCE   = 1000.0    # max gap (m) for two trace segments to be treated as connected

# --- Stage 1: surface meshing ---
DEPTH_CONTOUR_LEVELS = np.arange(0., 32000., 500.)  # depths (m) at which to draw contours
MESH_RESOLUTION      = 1000.0                         # target triangle size (m) for the raw surface
REBUILD_SURFACES     = True  # set True to regenerate a fault's surface even if its OBJ already exists
MERGE_UNCUT_MESHES   = True   # combine every Stage 1 VTK into MERGED_UNCUT_VTK for a quick inspection view

# --- Stage 2: cutting ---
MIN_CUT_DISTANCE = 2.0e3     # don't cut faults whose meshes are closer than this (m)
BOTTOM_DEPTH     = -31500.0  # depth (m) used when reasoning about cut extents
MERGE_CUT_MESHES = True      # combine every Stage 2 cut OBJ into MERGED_CUT_VTK for a quick inspection view

# --- Stage 2b: checks on the cut meshes ---
TRACE_CHECK_FLAG_DISTANCE = 2.0e3   # flag a cut mesh reaching this far (m) beyond the original trace
TRACE_CHECK_LATERAL_TOL   = 500.0   # surface trace within this (m) of the original counts as coincident
JAGGED_ANGLE_DEG          = 20.0    # direction change counting as a kink in the subsurface outline
JAGGED_MIN_SEGMENTS       = 3       # number of kinks needed to call a cut jagged
JAGGED_STEP_SCALES        = (500.0, 1000.0)  # chord lengths (m) at which to look for staircases
JAGGED_STEP_ANGLE         = 30.0    # direction change counting as a turn in a step
JAGGED_MIN_REVERSALS      = 3       # number of direction reversals needed to call a cut stepped

# --- Stage 3: MMG remeshing ---
TARGET_SIZE     = 1000.0     # uniform target edge length (m)
HAUSD           = 500.0      # max allowed deviation of the remesh from the original surface (m)
RIDGE_ANGLE_DEG = 45.0       # dihedral angle above which kinks are preserved

# --- Output directories ---
OBJ_DIR     = Path("test_objs")           # raw surfaces (OBJ) -- input to Stage 2
VTK_DIR     = Path("test_vtks")           # raw surfaces (VTK) -- for inspection
CONTOUR_DIR = Path("test_contours")       # depth contours (GeoJSON)
FINAL_DIR   = Path("test_final_meshes")   # cut surfaces -- input to Stage 3
REMESH_DIR  = Path("test_remeshed_vtks")  # remeshed surfaces
MERGED_VTK  = Path("merged_mesh_remeshed.vtk")
MERGED_UNCUT_VTK = Path("merged_mesh_uncut.vtk")  # combined Stage 1 raw surfaces, for inspection
MERGED_CUT_VTK   = Path("merged_mesh_cut.vtk")    # combined Stage 2 cut (unremeshed) surfaces, for inspection
TRACE_CHECK_CSV      = Path("surface_trace_check.csv")      # Stage 2b report
TRACE_CHECK_GEOJSON  = Path("surface_trace_check.geojson")  # Stage 2b geometry, for QGIS
JAGGED_CHECK_CSV     = Path("jagged_cut_check.csv")         # Stage 2b report
JAGGED_CHECK_GEOJSON = Path("jagged_cut_check.geojson")     # Stage 2b geometry, for QGIS
JAGGED_BLAME_CSV     = Path("jagged_cut_blame.csv")         # Stage 2b: which cut caused each kink
JAGGED_BLAME_GEOJSON = Path("jagged_cut_blame.geojson")     # Stage 2b geometry, for QGIS
CUTS_APPLIED_CSV     = Path("cuts_applied.csv")             # what was actually cut against what

SHOW_PLOT = False  # set True for an interactive 3-D view of the merged result at the end

for d in (OBJ_DIR, VTK_DIR, CONTOUR_DIR, FINAL_DIR, REMESH_DIR):
    d.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Build the fault network. This is the shared foundation for every stage:
# read the raw traces and reproject to NZTM, find connections between
# segments within DIST_TOLERANCE, read the curated fault systems (which
# segments form one continuous fault) and the cutting hierarchy (which fault
# truncates which). fault_data is built once and reused by all three stages.
# ---------------------------------------------------------------------------
fault_data = LeapfrogMultiFault.from_nz_cfm_shp(
    str(FAULT_SHP), remove_colons=True, epsg=EPSG,
    smoothing_n=None, dip_multiplier=1.0, exclude_zero=False)
fault_data.segment_distance_tolerance = DIST_TOLERANCE
fault_data.find_connections(verbose=False)

# Curated multi-segment faults, and the order in which faults cut each other.
fault_data.read_fault_systems(str(FAULT_SYSTEMS))
fault_data.generate_curated_faults()
fault_data.read_cutting_hierarchy(str(CUTTING_HIERARCHY))

print(f"{len(fault_data.curated_faults)} curated faults; "
      f"cutting hierarchy has {len(fault_data.cutting_hierarchy)} entries")

# ---------------------------------------------------------------------------
# Trace extensions (human-in-the-loop). Faults often need to be extended a
# little along strike so that neighbouring faults actually meet and cut
# cleanly. suggest_trace_extensions() writes proposed extensions for review
# in GIS (safe to re-run, never overwrites the edited file); once reviewed
# and saved as TRACE_EXT_EDITED, the curated result is applied to the shared
# fault_data before meshing, so the surfaces built in Stage 1 reflect them.
# No curated file exists yet for the whole-NZ dataset, so this step is
# skipped until TRACE_EXT_EDITED is created.
# ---------------------------------------------------------------------------
fault_data.suggest_trace_extensions(
    out_file="suggested_trace_extensions.csv",
    geojson_out_file="suggested_trace_extensions.geojson",
    fit_distance=5.e3, extend_distance=40.e3, proximity_threshold=1.e3)

if TRACE_EXT_EDITED.exists():
    fault_data.read_trace_extensions(str(TRACE_EXT_EDITED))
    # Extensions ruled out by review_trace_extensions.py are held in their own
    # file, so the curated list above is never rewritten by a tool.
    if TRACE_EXT_REMOVED.exists():
        fault_data.read_removed_trace_extensions(str(TRACE_EXT_REMOVED))
    fault_data.apply_trace_extensions()
else:
    print(f"Skipping trace extensions: {TRACE_EXT_EDITED} not found "
          f"(review suggested_trace_extensions.csv and save as {TRACE_EXT_EDITED} to enable).")

# ---------------------------------------------------------------------------
# Stage 1 - Build fault surfaces.
#
# For each curated fault we draw a stack of depth contours (every 500 m down
# to 32 km) then triangulate the surface spanning those contours. If the
# proper surface-meshing step fails for a fault (it can, for awkward
# geometries), fall back to a simpler contour-based mesh and save the
# offending contours to failed_contours_*.geojson for inspection, so one bad
# fault never stops the whole run. Each surface is written as OBJ (input to
# Stage 2), VTK (for quick viewing), and the contours as GeoJSON.
# ---------------------------------------------------------------------------
for fault in fault_data.curated_faults:
    obj_file = OBJ_DIR / f"{fault.name}_depth_contours.obj"
    if obj_file.exists() and not REBUILD_SURFACES:
        print(f"Skipping {fault.name}: surface already exists ({obj_file})")
        continue

    print(fault.name)
    try:
        fault.generate_depth_contours(DEPTH_CONTOUR_LEVELS, smoothing=False)
        mesh = fault.mesh_fault_surface(check_mesh=False, resolution=MESH_RESOLUTION)
    except Exception as e:
        # Fallback: simple contour mesh, and save the contours that tripped us up.
        print(f"  Failed to mesh {fault.name}: {e}")
        write_geojson_safe(fault.contours, f"failed_contours_{fault.name}.geojson", driver="GeoJSON")
        mesh = fault.mesh_simple_contours(DEPTH_CONTOUR_LEVELS)

    write_geojson_safe(fault.contours, CONTOUR_DIR / f"{fault.name}_depth_contours.geojson", driver="GeoJSON")
    mesh.write(str(VTK_DIR / f"{fault.name}_depth_contours.vtk"))
    mesh.write(str(OBJ_DIR / f"{fault.name}_depth_contours.obj"))

# Merge every raw Stage 1 surface into one VTK so the whole (uncut) fault
# network can be inspected in one go (e.g. in ParaView) before cutting.
if MERGE_UNCUT_MESHES:
    uncut_paths = sorted(VTK_DIR.glob("*.vtk"))
    uncut_meshes = [pv.from_meshio(meshio.read(str(p))) for p in uncut_paths]
    pv.merge(uncut_meshes).save(str(MERGED_UNCUT_VTK))
    print(f"Merged {len(uncut_meshes)} uncut surfaces -> {MERGED_UNCUT_VTK}")

# ---------------------------------------------------------------------------
# Stage 2 - Cut the surfaces.
#
# Truncate faults against each other and against the base depth surface, so
# the final meshes don't overlap or extend below the model. Load the depth
# raster as a PyVista surface (the lower bound everything is trimmed to) and
# the additional_cuts overrides (extra truncation relationships not captured
# by the hierarchy alone; an optional inverse, read_excluded_cuts(), exists
# for pairs that must never be cut -- both default to empty sets). Then read
# the Stage 1 OBJ surfaces back from disk into each fault's .mesh, which
# keeps this stage self-contained so you can re-cut without re-running
# Stage 1.
# ---------------------------------------------------------------------------
depth_pyvista = read_raster(str(DEPTH_RASTER), use_z=True, out_crs=f"EPSG:{EPSG}")
fault_data.read_additional_cuts(str(ADDITIONAL_CUTS))
# Override list of pairs that must NEVER be cut. Same 2-column, header-less CSV
# format as additional_cuts: each row is `fault_to_cut,cutting_fault`. Written
# by review_jagged_cuts.py, which walks the jagged cuts found below in 3-D and
# lets you rule out the ones that turn out to be wrong; skipped until that file
# exists, so the first run through needs nothing.
if EXCLUDED_CUTS.exists():
    fault_data.read_excluded_cuts(str(EXCLUDED_CUTS))
else:
    print(f"No {EXCLUDED_CUTS} yet: every cut is decided by the hierarchy and geometry.")

# Load the Stage 1 surfaces back onto each fault. FaultMesh.from_file names a
# mesh after its file, so rename each to its curated fault name: those names are
# what the cutting record below (and the CSVs) are keyed on.
for fault in fault_data.curated_faults:
    obj_file = OBJ_DIR / f"{fault.name}_depth_contours.obj"
    if obj_file.exists():
        fault.mesh = FaultMesh.from_file(str(obj_file))
        fault.mesh.name = fault.name

# Lookup of faults that actually have a mesh, used during cutting.
cutting_dict = {f.name: f for f in fault_data.curated_faults if f.mesh is not None}
print(f"{len(cutting_dict)} faults have meshes available for cutting")

# Apply the cutting hierarchy. Walk the faults in hierarchy order. For each
# fault the candidate cutters are every higher-priority fault, plus any
# fault named for it in additional_cuts -- so an additional cut overrides
# the hierarchy and applies even when the named cutter ranks lower (or isn't
# in the hierarchy at all). fancy_cutting=True uses the more robust
# intersection-following cut. Finally each fault is trimmed against the
# depth surface with cut_mesh_pv() and the result written as <name>_cut.obj
# -- the input to Stage 3.
threshold = 0.7  # currently unused by the library; kept for API compatibility

for fault_name in fault_data.cutting_hierarchy:
    print(f"Processing {fault_name} ...")
    fault = fault_data.name_dict[fault_name]
    if fault.mesh is None:
        continue

    # Faults ranked above this one are candidates to cut it.
    faults_to_cut = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(fault_name)]
    # additional_cuts override the hierarchy: a fault listed as cutting this one is
    # forced in as a candidate even when it sits lower in (or outside) the
    # hierarchy, where the loop above would otherwise never consider it.
    forced_cutters = [cut for (cut_fault, cut) in fault_data.additional_cuts
                      if cut_fault == fault_name and cut not in faults_to_cut]
    for cut_name in faults_to_cut + forced_cutters:
        if cut_name not in cutting_dict:
            continue
        cutting_mesh = cutting_dict[cut_name].mesh
        # Manual overrides win first. should_cut() matches on the curated fault
        # names and returns True/False to force/skip a cut, or None to fall
        # through to the geometric test below. additional_cuts/excluded_cuts
        # default to empty sets, so this is a no-op unless you loaded the CSVs.
        decision = fault_data.should_cut(fault_name, cut_name)
        if decision is None:
            # Only the geometric test needs the cutter's own higher-priority
            # context, so build it lazily (a forced cutter may sit outside the
            # hierarchy, where .index() would fail).
            higher = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(cut_name)]
            higher_meshes = [cutting_dict[n].mesh for n in higher if n in cutting_dict]
            decision = fault.mesh.decide_whether_to_cut(cutting_mesh, threshold=threshold,
                                                        min_distance=MIN_CUT_DISTANCE,
                                                        higher_meshes=higher_meshes,
                                                        bottom_depth=BOTTOM_DEPTH, fancy_cutting=True)
        if decision:
            print(f"  Cutting {fault.name} by {cut_name}")
            fragment = cutting_mesh.generate_cutting_mesh(fault.mesh, max_distance=5.e3)
            new_mesh = fault.mesh.cut_mesh(cutting_mesh,
                                           fault_trace=fault.original_nztm_trace_array,
                                           cutting_fragment=fragment, fancy_cutting=True)
            fault_data.name_dict[fault.name].mesh = new_mesh
            cutting_dict[fault.name].mesh = new_mesh

    # Trim against the base depth surface, then save the final cut mesh.
    try:
        depth_trimmed = fault_data.name_dict[fault.name].mesh.cut_mesh_pv(
            depth_pyvista, fault_trace=fault.original_nztm_trace_array)
        fault_data.name_dict[fault.name].mesh = depth_trimmed
    except Exception as e:
        print(f"  Error trimming depth for {fault.name}: {e}")
        continue
    fault_data.name_dict[fault.name].mesh.mesh.write(str(FINAL_DIR / f"{fault.name}_cut.obj"))
    print(f"  Done: {fault.name}")

# Record what was actually cut against what. Each mesh carries the list of
# surfaces that really changed it (FaultMesh.cut_against), which is a much
# smaller set than "everything ranked above it" -- Waitangi, for instance, is
# tested against 327 candidates and cut by exactly one. Stage 2b uses this to
# say for certain which cut left a jagged edge.
cut_records = [(fault.name, cutter)
               for fault in fault_data.curated_faults if fault.mesh is not None
               for cutter in fault.mesh.cut_against]
pd.DataFrame(cut_records, columns=["fault_name", "cutter"]).to_csv(CUTS_APPLIED_CSV, index=False)
print(f"Recorded {len(cut_records)} applied cuts -> {CUTS_APPLIED_CSV}")

# Merge every cut (but not yet remeshed) surface into one VTK so the whole
# fault network can be inspected in one go (e.g. in ParaView) before Stage 3.
if MERGE_CUT_MESHES:
    cut_paths = sorted(FINAL_DIR.glob("*_cut.obj"))
    cut_meshes = [pv.from_meshio(meshio.read(str(p))) for p in cut_paths]
    pv.merge(cut_meshes).save(str(MERGED_CUT_VTK))
    print(f"Merged {len(cut_meshes)} cut surfaces -> {MERGED_CUT_VTK}")

# ---------------------------------------------------------------------------
# Stage 2b - Check the cut meshes. Both checks read the cut meshes back from
# disk, so they can be re-run on their own, and both write a CSV plus a GeoJSON
# holding only the offending geometry, to be loaded over the traces in QGIS.
#
# 1. Surface trace vs original trace. Trace extensions are there to grow a
#    fault's *subsurface* area so it reaches its neighbours and can be cut
#    cleanly -- they are not meant to survive at the surface. This takes each
#    mesh's surface trace (the outline at z = 0) and measures how far it runs
#    past the end of the original, un-extended trace. Anything over
#    TRACE_CHECK_FLAG_DISTANCE is flagged, to be fixed either by forcing the
#    missing truncation in additional_cuts.csv or by dropping the extension
#    from edited_trace_extensions.csv.
#
# 2. Jagged cuts. A cut that has gone wrong leaves a serrated outline, zig-
#    zagging between the triangles it kept and the ones it dropped. This walks
#    the subsurface outline (the surface trace is excluded: a real trace is
#    expected to bend sharply) and counts segments whose direction differs from
#    both their neighbours by JAGGED_ANGLE_DEG or more. Each kink is then
#    blamed on whichever cut left it, by finding the nearest surface that was
#    entitled to cut that fault -- so a jagged fault comes with the name of the
#    cut to revisit rather than just a score.
# ---------------------------------------------------------------------------
trace_check = fault_data.check_surface_trace_extensions(
    str(FINAL_DIR), flag_distance=TRACE_CHECK_FLAG_DISTANCE,
    lateral_tolerance=TRACE_CHECK_LATERAL_TOL,
    out_file=str(TRACE_CHECK_CSV), geojson_out_file=str(TRACE_CHECK_GEOJSON))

jagged_check = fault_data.check_jagged_cuts(
    str(FINAL_DIR), angle_threshold=JAGGED_ANGLE_DEG,
    min_kinked_segments=JAGGED_MIN_SEGMENTS, step_scales=JAGGED_STEP_SCALES,
    step_angle=JAGGED_STEP_ANGLE, min_reversals=JAGGED_MIN_REVERSALS,
    out_file=str(JAGGED_CHECK_CSV), geojson_out_file=str(JAGGED_CHECK_GEOJSON))

# cuts_applied.csv turns the blame from a guess into a fact; without it (Stage 2
# not run in this session) attribution falls back to re-testing each candidate.
jagged_blame = fault_data.attribute_jagged_cuts(
    jagged_check[jagged_check["flagged"]], str(OBJ_DIR), depth_surface=depth_pyvista,
    actual_cuts=str(CUTS_APPLIED_CSV) if CUTS_APPLIED_CSV.exists() else None,
    min_cut_distance=MIN_CUT_DISTANCE, bottom_depth=BOTTOM_DEPTH, cut_threshold=threshold,
    out_file=str(JAGGED_BLAME_CSV), geojson_out_file=str(JAGGED_BLAME_GEOJSON))

# # ---------------------------------------------------------------------------
# # Stage 3 - Remesh with MMG.
# #
# # The cut surfaces have uneven, often sliver-shaped triangles left over from
# # the contouring and cutting. Clean each one into uniform, near-equilateral
# # triangles of TARGET_SIZE using MMG (mmgs). MMG remeshes the surface in
# # place (local edge split/collapse/swap, nodes moved back onto the original
# # surface) rather than flattening it onto a plane and re-triangulating --
# # for a long, warped fault sheet like the Alpine Fault that flattening folds
# # over itself and silently drops the folded region (about half of the
# # Alpine Fault was being deleted by a parametrise-and-regenerate remesher).
# # HAUSD controls how tightly the new mesh hugs the original surface;
# # RIDGE_ANGLE_DEG tells MMG which sharp kinks to preserve rather than smooth.
# # ---------------------------------------------------------------------------
# remeshed_paths = []
# for obj_path in sorted(FINAL_DIR.glob("*_cut.obj")):
#     name = obj_path.stem.removesuffix("_cut")
#     out_path = REMESH_DIR / f"{name}_remeshed.vtk"
#     print(f"Remeshing {name}")
#     try:
#         fault_mesh = FaultMesh.from_file(str(obj_path))
#         remeshed = fault_mesh.remesh(target_size=TARGET_SIZE, hausd=HAUSD,
#                                      ridge_angle_deg=RIDGE_ANGLE_DEG, verbose=False)
#         meshio.write(str(out_path), remeshed.mesh)
#         remeshed_paths.append(out_path)
#     except Exception as e:
#         print(f"  failed: {e}")

# # Merge every remeshed surface into one VTK so the whole fault network can be
# # opened and checked in one go (e.g. in ParaView).
# meshes = [pv.from_meshio(meshio.read(str(p))) for p in remeshed_paths]
# combined = pv.merge(meshes)
# combined.save(str(MERGED_VTK))
# print(f"Merged {len(meshes)} remeshed surfaces -> {MERGED_VTK}")

# # (Optional) Quick 3-D view of the merged result.
# if SHOW_PLOT:
#     combined = pv.read(str(MERGED_VTK))
#     plotter = pv.Plotter()
#     plotter.add_mesh(combined, show_edges=True, color="lightgray")
#     plotter.show()
