"""Interactive review of the jagged cuts flagged by Stage 2b of the pipeline.

Works through the faults in jagged_cut_blame.csv, worst first. For each one it
opens a 3-D PyVista view of the jagged fault together with every fault that may
have cut it -- each a semi-transparent surface with its triangle edges in its
own colour -- and highlights the rough parts of the outline. Close the window
and the terminal offers the decision: leave it alone, or write the suspected cut
to EXCLUDED_CUTS so the next pipeline run does not make it.

Run from this file's directory (scratch/whole_nz_no_leapfrog/) inside the
leapfrog-fault-models conda environment, after Stage 2b has produced
jagged_cut_check.csv, jagged_cut_blame.csv and jagged_cut_blame.geojson:

    python review_jagged_cuts.py                 # review everything not yet done
    python review_jagged_cuts.py --only Waitangi # one fault
    python review_jagged_cuts.py --all           # start again from the top

Decisions are appended to EXCLUDED_CUTS and REVIEW_LOG as they are made, so the
review can be stopped at any point (q, or Ctrl-C) and picked up later.
"""
import argparse
import datetime
import sys
from pathlib import Path

import pandas as pd
import geopandas as gpd

from review_common import (FaultView, FAULT_COLOUR, append_row, load_logged, load_pairs,
                           mesh_path)

# --- Inputs (relative to this script's folder) ---
JAGGED_CHECK_CSV = Path("jagged_cut_check.csv")          # roughness scores per fault
JAGGED_BLAME_CSV = Path("jagged_cut_blame.csv")          # which cut each rough part belongs to
JAGGED_BLAME_GEOJSON = Path("jagged_cut_blame.geojson")  # geometry of the rough parts

# --- Outputs ---
EXCLUDED_CUTS = Path("excluded_cuts.csv")    # pairs that must never be cut
REVIEW_LOG = Path("jagged_review_log.csv")   # what has been reviewed, and when

NOT_A_FAULT = ("unattributed", "base depth surface")  # cannot be excluded


def show_fault(fault_name, attributions, kink_geometries, z_scale, prefer_cut, screenshot=None):
    """Open the 3-D view of one jagged fault and the faults that may have cut it."""
    view = FaultView(fault_name, screenshot=screenshot)

    path = mesh_path(fault_name, prefer_cut=True)
    if path is None:
        print(f"  No mesh found for {fault_name}")
        return
    view.add_surface(path, f"{fault_name} (jagged)", colour=FAULT_COLOUR, focus=True)

    for _, row in attributions.iterrows():
        if row["cutter"] in NOT_A_FAULT:
            continue
        cutter_path = mesh_path(row["cutter"], prefer_cut=prefer_cut)
        if cutter_path is None:
            print(f"  No mesh found for {row['cutter']}")
            continue
        view.add_surface(cutter_path, row["cutter"])

    view.add_lines(kink_geometries, "rough outline")
    if z_scale != 1.0:
        view.plotter.set_scale(zscale=z_scale)
    view.show()


def describe(fault_name, attributions, check_row, already_excluded):
    """Print the attribution table and return the list of excludable cutters."""
    print(f"\n{'=' * 78}")
    if check_row is not None:
        summary = []
        if check_row.get("n_reversals", 0):
            summary.append(f"{int(check_row['n_reversals'])} step reversals at "
                           f"{check_row['reversal_scale_m']:.0f} m "
                           f"(steps of {check_row['step_amplitude_m']:.0f} m)")
        if check_row.get("n_kinked_segments", 0):
            summary.append(f"{int(check_row['n_kinked_segments'])} kinked segments "
                           f"(worst turn {check_row['max_turn_deg']:.0f} deg)")
        print(f"{fault_name}: {'; '.join(summary)}")
    else:
        print(fault_name)

    options = []
    for _, row in attributions.iterrows():
        cutter = row["cutter"]
        confirmed = row.get("cut_confirmed", "")
        if pd.isna(confirmed):
            confirmed = ""
        detail = (f"{row['n_segments']:3d} rough segments, "
                  f"median {row['median_distance_m']:6.0f} m"
                  f"{'' if confirmed == '' else f', cut confirmed: {confirmed}'}")
        if cutter in NOT_A_FAULT:
            print(f"     -  {cutter:<38s} {detail}")
            continue
        options.append(cutter)
        flag = "  [already excluded]" if (fault_name, cutter) in already_excluded else ""
        print(f"    {len(options)})  {cutter:<38s} {detail}{flag}")
    return options


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", metavar="FAULT", help="review just this fault")
    parser.add_argument("--start-from", metavar="FAULT", help="skip ahead to this fault")
    parser.add_argument("--all", action="store_true",
                        help="include faults already in the review log")
    parser.add_argument("--z-scale", type=float, default=1.0,
                        help="vertical exaggeration of the 3-D view (default 1)")
    parser.add_argument("--uncut-cutters", action="store_true",
                        help="show the cutting faults' raw surfaces rather than their cut ones")
    parser.add_argument("--screenshots", metavar="DIR",
                        help="render every fault to a PNG in DIR instead of reviewing them")
    args = parser.parse_args()

    for required in (JAGGED_BLAME_CSV, JAGGED_BLAME_GEOJSON):
        if not required.exists():
            sys.exit(f"{required} not found -- run Stage 2b of the pipeline first.")

    blame = pd.read_csv(JAGGED_BLAME_CSV)
    kinks = gpd.read_file(JAGGED_BLAME_GEOJSON)
    check = pd.read_csv(JAGGED_CHECK_CSV).set_index("fault_name") \
        if JAGGED_CHECK_CSV.exists() else None

    # Keep the order of the blame table, which is worst-first.
    fault_names = list(dict.fromkeys(blame["fault_name"]))
    if args.only:
        if args.only not in fault_names:
            sys.exit(f"{args.only} is not in {JAGGED_BLAME_CSV}")
        fault_names = [args.only]
    else:
        if args.start_from:
            if args.start_from not in fault_names:
                sys.exit(f"{args.start_from} is not in {JAGGED_BLAME_CSV}")
            fault_names = fault_names[fault_names.index(args.start_from):]
        if not args.all:
            reviewed = load_logged(REVIEW_LOG)
            skipped = [name for name in fault_names if name in reviewed]
            fault_names = [name for name in fault_names if name not in reviewed]
            if skipped:
                print(f"Skipping {len(skipped)} faults already in {REVIEW_LOG} "
                      f"(use --all to include them)")

    if not fault_names:
        print("Nothing left to review.")
        return

    if args.screenshots:
        directory = Path(args.screenshots)
        directory.mkdir(parents=True, exist_ok=True)
        for fault_name in fault_names:
            out_path = directory / f"{fault_name}.png"
            show_fault(fault_name, blame[blame["fault_name"] == fault_name],
                       kinks[kinks["fault_name"] == fault_name]["geometry"],
                       args.z_scale, prefer_cut=not args.uncut_cutters,
                       screenshot=str(out_path))
            print(f"  {out_path}")
        print(f"\nRendered {len(fault_names)} faults to {directory}")
        return

    print(f"{len(fault_names)} faults to review. "
          f"Excluded cuts are appended to {EXCLUDED_CUTS} as you go.")

    for position, fault_name in enumerate(fault_names, start=1):
        attributions = blame[blame["fault_name"] == fault_name]
        fault_kinks = kinks[kinks["fault_name"] == fault_name]["geometry"]
        check_row = check.loc[fault_name] if check is not None and fault_name in check.index \
            else None

        print(f"\n[{position}/{len(fault_names)}]", end="")
        options = describe(fault_name, attributions, check_row, load_pairs(EXCLUDED_CUTS))

        while True:
            show_fault(fault_name, attributions, fault_kinks, args.z_scale,
                       prefer_cut=not args.uncut_cutters)
            prompt = ("\n  [Enter] leave alone"
                      f"{'   [1-%d] exclude that cut' % len(options) if options else ''}"
                      "   [r] show again   [q] quit\n  > ")
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
                append_row(REVIEW_LOG, [fault_name, "no change", "",
                                        datetime.datetime.now().isoformat(timespec="seconds")],
                           header=["fault_name", "decision", "excluded_cutter", "reviewed_at"])
                print("  left alone")
                break

            # One or more numbers, e.g. "2" or "1,3".
            try:
                chosen = [options[int(part) - 1] for part in answer.split(",")]
            except (ValueError, IndexError):
                print(f"  Did not understand '{answer}'")
                continue

            already_excluded = load_pairs(EXCLUDED_CUTS)
            for cutter in chosen:
                if (fault_name, cutter) in already_excluded:
                    print(f"  {fault_name} / {cutter} is already in {EXCLUDED_CUTS}")
                    continue
                append_row(EXCLUDED_CUTS, [fault_name, cutter])
                append_row(REVIEW_LOG,
                           [fault_name, "excluded cut", cutter,
                            datetime.datetime.now().isoformat(timespec="seconds")],
                           header=["fault_name", "decision", "excluded_cutter", "reviewed_at"])
                print(f"  excluded: {fault_name} is no longer cut by {cutter}")
            break

    print(f"\nReviewed {len(fault_names)} faults. "
          f"Re-run the pipeline to rebuild the meshes with {EXCLUDED_CUTS} applied.")


if __name__ == "__main__":
    main()
