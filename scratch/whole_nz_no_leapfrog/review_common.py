"""Shared pieces of the Stage 2b review tools.

review_jagged_cuts.py and review_trace_extensions.py both show a fault in 3-D
next to the faults around it and then ask what to do about it, so the viewer,
the colour scheme and the decision-file plumbing live here.
"""
import csv
from pathlib import Path

import numpy as np
import pyvista as pv
from matplotlib.colors import to_rgb

# --- Mesh directories (relative to the pipeline folder) ---
FINAL_DIR = Path("test_final_meshes")   # cut surfaces
OBJ_DIR = Path("test_objs")             # uncut surfaces

# --- Display ---
# The fault under review always takes the first colour; the faults around it
# take the rest, so each surface and its triangle edges are distinguishable.
FAULT_COLOUR = "#d1495b"
OTHER_COLOURS = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd",
                 "#17becf", "#8c564b", "#7f7f7f", "#bcbd22"]
HIGHLIGHT_COLOUR = "#ffd400"
SURFACE_OPACITY = 0.45
EDGE_WIDTH = 1
EDGE_OPACITY = 0.35
# The camera is framed on the fault under review, expanded by this much, so a
# neighbour several times its length does not shrink it to a speck.
FRAME_MARGIN = 0.35
WINDOW_SIZE = (1400, 900)
# The legend box scales its text to its own height, so the height is set from
# the number of entries and long fault names are shortened to fit the width.
LEGEND_ROW_HEIGHT = 0.022
LEGEND_LABEL_CHARS = 20


def darker(colour, factor=0.55):
    """A darker shade of `colour`, for drawing that surface's triangle edges."""
    return tuple(channel * factor for channel in to_rgb(colour))


def mesh_path(name, prefer_cut=True):
    """Locate a fault's surface, preferring the cut version over the raw one."""
    candidates = [FINAL_DIR / f"{name}_cut.obj", OBJ_DIR / f"{name}_depth_contours.obj"]
    if not prefer_cut:
        candidates.reverse()
    for path in candidates:
        if path.exists():
            return path
    return None


def scene_span(bounds):
    """Diagonal length of a (xmin, xmax, ymin, ymax, zmin, zmax) box."""
    bounds = np.array(bounds, dtype=float)
    return max(float(np.linalg.norm(bounds[1::2] - bounds[0::2])), 1.0)


def lines_to_polydata(geometries):
    """Build a PyVista line set from shapely (Multi)LineStrings."""
    points, lines = [], []
    for geometry in geometries:
        parts = list(geometry.geoms) if hasattr(geometry, "geoms") else [geometry]
        for part in parts:
            coords = np.asarray(part.coords)
            if coords.shape[1] == 2:  # surface traces come through flat
                coords = np.column_stack([coords, np.zeros(len(coords))])
            start = len(points)
            points.extend(coords)
            lines.extend([len(coords)] + list(range(start, start + len(coords))))
    if not points:
        return None
    return pv.PolyData(np.array(points), lines=np.array(lines))


class FaultView:
    """A 3-D view of one fault and its neighbours, built up a surface at a time."""

    def __init__(self, title, screenshot=None):
        self.plotter = pv.Plotter(window_size=WINDOW_SIZE, off_screen=screenshot is not None)
        self.plotter.add_text(title if screenshot else
                              f"{title} -- close this window to choose what to do",
                              font_size=11)
        self.screenshot = screenshot
        self.n_labels = 0
        self.focus_surface = None
        self._colour_index = 0

    def add_surface(self, path, label, colour=None, focus=False):
        """Draw one fault as a translucent surface with its triangles outlined.

        The outlines go on as their own actor: edges drawn as part of a
        translucent surface come out translucent too, which at these scales
        leaves the triangulation invisible.
        """
        if colour is None:
            colour = OTHER_COLOURS[self._colour_index % len(OTHER_COLOURS)]
            self._colour_index += 1
        surface = pv.read(str(path))
        self.plotter.add_mesh(surface, color=colour, opacity=SURFACE_OPACITY,
                              label=label[:LEGEND_LABEL_CHARS])
        self.plotter.add_mesh(surface, style="wireframe", color=darker(colour),
                              line_width=EDGE_WIDTH, opacity=EDGE_OPACITY)
        self.n_labels += 1
        if focus:
            self.focus_surface = surface
        return surface

    def add_lines(self, geometries, label, colour=HIGHLIGHT_COLOUR, width=8):
        """Draw a set of lines as thick tubes, e.g. the part of an outline at fault."""
        poly = lines_to_polydata(geometries)
        if poly is None:
            return
        self.plotter.add_mesh(poly, color=colour, line_width=width,
                              render_lines_as_tubes=True, label=label)
        self.n_labels += 1

    def show(self):
        """Frame the view on the fault under review and display it."""
        self.plotter.add_legend(bcolor="white", loc="upper right",
                                size=(0.24, max(2, self.n_labels) * LEGEND_ROW_HEIGHT))
        self.plotter.add_axes()
        self.plotter.view_isometric()

        # Frame on the fault being reviewed rather than on everything drawn: its
        # neighbours are often several times longer and would fill the view. The
        # camera is moved rather than reset, since a reset gets undone by the
        # render that follows; sliding focal point and position together keeps
        # the isometric direction and just re-centres.
        if self.focus_surface is not None:
            focus = np.array(self.focus_surface.center, dtype=float)
            shift = focus - np.array(self.plotter.camera.focal_point, dtype=float)
            self.plotter.camera.focal_point = focus
            self.plotter.camera.position = (np.array(self.plotter.camera.position, dtype=float)
                                            + shift)
            self.plotter.camera.zoom(scene_span(self.plotter.bounds) /
                                     (scene_span(self.focus_surface.bounds)
                                      * (1. + 2. * FRAME_MARGIN)))
        self.plotter.show(screenshot=self.screenshot)


def append_row(path, row, header=None):
    """Append one row to a CSV, creating it (with `header`) if it does not exist."""
    new_file = not Path(path).exists()
    with open(path, "a", newline="") as handle:
        writer = csv.writer(handle)
        if new_file and header is not None:
            writer.writerow(header)
        writer.writerow(row)


def load_pairs(path):
    """Read a header-less two-column CSV as a set of tuples."""
    import pandas as pd
    if not Path(path).exists():
        return set()
    frame = pd.read_csv(path, header=None, names=["a", "b"])
    return {(str(a).strip(), str(b).strip()) for a, b in zip(frame["a"], frame["b"])}


def load_logged(path, column="fault_name"):
    """Values already recorded in a review log, for resuming a session."""
    import pandas as pd
    if not Path(path).exists():
        return set()
    return set(pd.read_csv(path)[column])
