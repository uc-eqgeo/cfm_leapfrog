"""Interactive review of the trace-extension overshoots found by Stage 2b.

A trace extension exists to grow a fault's *subsurface* area so it reaches its
neighbours and can be cut cleanly. It is not meant to survive at the surface.
Where surface_trace_check.csv says a cut mesh still reaches the surface well
past the end of the original trace, one of two things has gone wrong:

  * the fault was expected to terminate against a neighbour and does not, so the
    fix is an additional cut; or
  * the extension was not needed in the first place, so the fix is to drop it.

For each flagged end this opens a 3-D view of the fault, the overshooting piece
of trace, and the faults it was expected to terminate against, then asks which
of those two it is.

Run from this file's directory (scratch/whole_nz_no_leapfrog/) inside the
leapfrog-fault-models conda environment, after Stage 2b has produced
surface_trace_check.csv and surface_trace_check.geojson:

    python review_trace_extensions.py                    # everything not yet done
    python review_trace_extensions.py --only Awakeri     # one fault
    python review_trace_extensions.py --all              # start again from the top

Decisions append to ADDITIONAL_CUTS and REMOVED_EXTENSIONS as they are made, so
the review can be stopped at any point (q, or Ctrl-C) and picked up later. The
curated edited_trace_extensions.csv is never modified: removals go to their own
file, which the pipeline subtracts after reading the curated one.
"""
import argparse
import datetime
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree

from review_common import (FaultView, FAULT_COLOUR, HIGHLIGHT_COLOUR, append_row,
                           load_logged, load_pairs, mesh_path)

# --- Inputs (relative to this script's folder) ---
TRACE_CHECK_CSV = Path("surface_trace_check.csv")          # overshoot per fault end
TRACE_CHECK_GEOJSON = Path("surface_trace_check.geojson")  # geometry of the overshoots
FAULT_SHP = Path("NZ_CFM_v1_0_rs1km_modified_GDM.gpkg")
FAULT_SYSTEMS = Path("NZ_CFM_v1_0_rs1km_modified_GDM_connected_edited.csv")
CUTTING_HIERARCHY = Path("NZ_CFM_v1_0_rs1km_modified_GDM_hierarchy_edited.csv")

# --- Outputs ---
ADDITIONAL_CUTS = Path("additional_cuts.csv")               # cuts to force
REMOVED_EXTENSIONS = Path("removed_trace_extensions.csv")   # extensions to drop
REVIEW_LOG = Path("trace_review_log.csv")                   # what has been reviewed

EPSG = 2193
DIST_TOLERANCE = 1000.0
CANDIDATE_DISTANCE = 5.0e3   # how close a fault must pass to the overshoot to be a candidate
MAX_CANDIDATES = 6


def build_fault_data():
    """The fault network, needed for the connection graph behind find_terminations()."""
    from fault_mesh.faults.leapfrog import LeapfrogMultiFault

    fault_data = LeapfrogMultiFault.from_nz_cfm_shp(
        str(FAULT_SHP), remove_colons=True, epsg=EPSG,
        smoothing_n=None, dip_multiplier=1.0, exclude_zero=False)
    fault_data.segment_distance_tolerance = DIST_TOLERANCE
    fault_data.find_connections(verbose=False)
    fault_data.read_fault_systems(str(FAULT_SYSTEMS))
    fault_data.generate_curated_faults()
    fault_data.read_cutting_hierarchy(str(CUTTING_HIERARCHY))
    return fault_data


def overshoot_points(geometry):
    """The overshooting trace as an (n, 2) array of map coordinates."""
    parts = list(geometry.geoms) if hasattr(geometry, "geoms") else [geometry]
    return np.vstack([np.asarray(part.coords)[:, :2] for part in parts])


def find_candidates(fault_data, fault_name, overshoot):
    """Faults this one may have been meant to terminate against, nearest first.

    Two sources, unioned. find_terminations() gives the faults the connection
    graph says it runs into and that outrank it -- exactly "expected to
    terminate against". Geometry then adds anything whose trace passes close to
    the overshooting piece, which catches the cases the graph missed (and those
    are the interesting ones, since a missing connection is one reason the cut
    never happened).
    """
    from itertools import chain as ichain

    expected = set()
    try:
        expected = {name for name in ichain(*fault_data.find_terminations(fault_name))
                    if name != fault_name}
    except ValueError:
        pass  # not in the cutting hierarchy

    tree = cKDTree(overshoot)
    candidates = {}
    for other in fault_data.curated_faults:
        if other.name == fault_name:
            continue
        coords = np.asarray(other.nztm_trace.coords)[:, :2]
        distance = float(tree.query(coords)[0].min())
        if distance <= CANDIDATE_DISTANCE or other.name in expected:
            candidates[other.name] = (distance, other.name in expected)
    return sorted(candidates.items(), key=lambda item: item[1][0])[:MAX_CANDIDATES]


def show_end(fault_name, overshoot_geometries, candidates, original_trace, z_scale,
             screenshot=None):
    """3-D view of the fault, its overshooting trace, and the candidate terminators."""
    view = FaultView(fault_name, screenshot=screenshot)

    path = mesh_path(fault_name, prefer_cut=True)
    if path is None:
        print(f"  No mesh found for {fault_name}")
        return
    view.add_surface(path, f"{fault_name} (overshoots)", colour=FAULT_COLOUR, focus=True)

    for name, _ in candidates:
        candidate_path = mesh_path(name, prefer_cut=True)
        if candidate_path is None:
            print(f"  No mesh found for {name}")
            continue
        view.add_surface(candidate_path, name)

    if original_trace is not None:
        view.add_lines([original_trace], "original trace", colour="#333333", width=4)
    view.add_lines(overshoot_geometries, "extra surface trace", colour=HIGHLIGHT_COLOUR)

    if z_scale != 1.0:
        view.plotter.set_scale(zscale=z_scale)
    view.show()


def describe(row, candidates, already_cut, already_removed):
    """Print the overshoot and its candidate terminators; return the candidate names."""
    print(f"\n{'=' * 78}")
    print(f"{row['fault_name']} ({row['end']} end): "
          f"{row['max_extra_m'] / 1e3:.1f} km beyond the original trace, "
          f"{row['extra_length_m'] / 1e3:.1f} km of extra trace")
    if row.get("has_extension"):
        print(f"    a {row['extend_distance'] / 1e3:.0f} km trace extension is applied at this end")
    else:
        print("    no trace extension is recorded at this end -- check the cutting hierarchy")

    names = []
    for name, (distance, expected) in candidates:
        names.append(name)
        note = "expected termination" if expected else "nearby"
        flag = "  [cut already forced]" if (row["fault_name"], name) in already_cut else ""
        print(f"    {len(names)})  {name:<38s} {distance / 1e3:5.1f} km away, {note}{flag}")
    if not names:
        print("    (no candidate faults near this end)")

    if (row["fault_name"], str(row["end"])) in already_removed:
        print("    [extension already marked for removal]")
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", metavar="FAULT", help="review just this fault")
    parser.add_argument("--start-from", metavar="FAULT", help="skip ahead to this fault")
    parser.add_argument("--all", action="store_true",
                        help="include ends already in the review log")
    parser.add_argument("--include-unflagged", action="store_true",
                        help="also review overshoots below the flagging threshold")
    parser.add_argument("--z-scale", type=float, default=1.0,
                        help="vertical exaggeration of the 3-D view (default 1)")
    parser.add_argument("--screenshots", metavar="DIR",
                        help="render every end to a PNG in DIR instead of reviewing them")
    args = parser.parse_args()

    for required in (TRACE_CHECK_CSV, TRACE_CHECK_GEOJSON):
        if not required.exists():
            sys.exit(f"{required} not found -- run Stage 2b of the pipeline first.")

    checks = pd.read_csv(TRACE_CHECK_CSV)
    if not args.include_unflagged:
        checks = checks[checks["flagged"]]
    geometries = gpd.read_file(TRACE_CHECK_GEOJSON)

    if args.only:
        checks = checks[checks["fault_name"] == args.only]
        if not len(checks):
            sys.exit(f"{args.only} has no flagged overshoot in {TRACE_CHECK_CSV}")
    elif args.start_from:
        names = list(checks["fault_name"])
        if args.start_from not in names:
            sys.exit(f"{args.start_from} is not in {TRACE_CHECK_CSV}")
        checks = checks.iloc[names.index(args.start_from):]

    # An end is identified by fault and compass direction, so both go in the log key.
    if not args.all and not args.only:
        reviewed = load_logged(REVIEW_LOG, column="key")
        keep = [f"{name}|{end}" not in reviewed
                for name, end in zip(checks["fault_name"], checks["end"])]
        skipped = len(keep) - sum(keep)
        checks = checks[keep]
        if skipped:
            print(f"Skipping {skipped} ends already in {REVIEW_LOG} (use --all to include them)")

    if not len(checks):
        print("Nothing left to review.")
        return

    print(f"Building the fault network for the connection graph ...")
    fault_data = build_fault_data()

    if args.screenshots:
        directory = Path(args.screenshots)
        directory.mkdir(parents=True, exist_ok=True)

    print(f"{len(checks)} ends to review. Decisions are appended to {ADDITIONAL_CUTS} "
          f"and {REMOVED_EXTENSIONS} as you go.")

    for position, (_, row) in enumerate(checks.iterrows(), start=1):
        fault_name, end = row["fault_name"], str(row["end"])
        match = geometries[(geometries["fault_name"] == fault_name)
                           & (geometries["end"].astype(str) == end)]["geometry"]
        if not len(match):
            print(f"  No overshoot geometry for {fault_name} ({end}); skipping")
            continue

        candidates = find_candidates(fault_data, fault_name, overshoot_points(match.iloc[0]))
        fault = fault_data.name_dict.get(fault_name)
        original = fault.original_nztm_trace if fault is not None else None

        if args.screenshots:
            out_path = Path(args.screenshots) / f"{fault_name}_{end}.png"
            show_end(fault_name, match, candidates, original, args.z_scale,
                     screenshot=str(out_path))
            print(f"  {out_path}")
            continue

        print(f"\n[{position}/{len(checks)}]", end="")
        names = describe(row, candidates, load_pairs(ADDITIONAL_CUTS),
                         load_logged_removals())

        while True:
            show_end(fault_name, match, candidates, original, args.z_scale)
            prompt = ("\n  [Enter] leave alone"
                      f"{'   [1-%d] add cut against that fault' % len(names) if names else ''}"
                      "   [x] remove the trace extension   [r] show again   [q] quit\n  > ")
            try:
                answer = input(prompt).strip().lower()
            except (EOFError, KeyboardInterrupt):
                print("\nStopped.")
                return

            if answer == "q":
                print("Stopped.")
                return
            if answer == "r":
                continue

            if answer == "":
                log(fault_name, end, "no change", "")
                print("  left alone")
                break

            if answer == "x":
                if (fault_name, end) in load_logged_removals():
                    print(f"  {fault_name} ({end}) is already in {REMOVED_EXTENSIONS}")
                else:
                    append_row(REMOVED_EXTENSIONS, [fault_name, end],
                               header=["fault_name", "end"])
                    log(fault_name, end, "removed extension", "")
                    print(f"  the {end} trace extension on {fault_name} will not be applied")
                break

            try:
                chosen = [names[int(part) - 1] for part in answer.split(",")]
            except (ValueError, IndexError):
                print(f"  Did not understand '{answer}'")
                continue

            already_cut = load_pairs(ADDITIONAL_CUTS)
            for cutter in chosen:
                if (fault_name, cutter) in already_cut:
                    print(f"  {fault_name} / {cutter} is already in {ADDITIONAL_CUTS}")
                    continue
                append_row(ADDITIONAL_CUTS, [fault_name, cutter])
                log(fault_name, end, "added cut", cutter)
                print(f"  added: {fault_name} will be cut by {cutter}")
            break

    print(f"\nReviewed {len(checks)} ends. Re-run the pipeline to rebuild the meshes.")


def load_logged_removals():
    """Extension removals already recorded, as (fault_name, end) tuples."""
    if not REMOVED_EXTENSIONS.exists():
        return set()
    frame = pd.read_csv(REMOVED_EXTENSIONS)
    return {(str(name).strip(), str(end).strip())
            for name, end in zip(frame["fault_name"], frame["end"])}


def log(fault_name, end, decision, detail):
    """Record a decision so the review can be resumed."""
    append_row(REVIEW_LOG,
               [f"{fault_name}|{end}", fault_name, end, decision, detail,
                datetime.datetime.now().isoformat(timespec="seconds")],
               header=["key", "fault_name", "end", "decision", "detail", "reviewed_at"])


if __name__ == "__main__":
    main()
