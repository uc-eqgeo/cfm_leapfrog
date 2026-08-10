"""
Classes that implement the Leapfrog fault model. Inherit from GenericFault and GenericMultiFault.
"""
from __future__ import annotations
import os
from typing import List, Union
from itertools import product, chain

import xml.etree.ElementTree as ElemTree
from xml.dom import minidom

import numpy as np
import geopandas as gpd
import pandas as pd
import meshio
import pyvista as pv
from scipy.spatial import cKDTree

from shapely.affinity import translate
from shapely.ops import unary_union, linemerge
from shapely.geometry import LineString, MultiLineString, Point, Polygon, MultiPoint

from fault_mesh.faults.generic import GenericMultiFault, GenericFault, normalize_bearing, smallest_difference, calculate_strike_direction
from fault_mesh.utilities.smoothing import smooth_trace
from fault_mesh.utilities.cutting import cut_line_between_two_points, cut_line_at_point
from fault_mesh.utilities.graph import connected_nodes, suggest_combined_name
from fault_mesh.faults.connected import ConnectedFaultSystem
from fault_mesh.utilities.merging import merge_multiple_nearly_adjacent_segments

from fault_mesh.utilities.meshing import triangulate_contours
from fault_mesh.utilities.splines import spline_fit_contours
from fault_mesh.utilities.opensha import fault_trace_xml, fault_polygon_xml
from fault_mesh.faults.mesh import FaultMesh


class LeapfrogMultiFault(GenericMultiFault):
    def __init__(self, fault_geodataframe: gpd.GeoDataFrame, sort_sr: bool = False,
                 segment_distance_tolerance: float = 100., smoothing_n: int = None,
                 remove_colons: bool = True, dip_choice: str = "pref",
                 trimming_gradient: float = 1., epsg: int = None, dip_multiplier: float = 1.,
                 strike_multiplier: float = 0.5, check_optional_fields: bool = True):

        self._smoothing_n = smoothing_n
        super(LeapfrogMultiFault, self).__init__(fault_geodataframe=fault_geodataframe, 
                                                 sort_sr=sort_sr,
                                                 remove_colons=remove_colons, dip_choice=dip_choice,
                                                 check_optional_fields=check_optional_fields)

        self._segment_distance_tolerance = segment_distance_tolerance
        self._cutting_hierarchy = []
        self._hierarchy_dict = None
        self._includes_connected = False
        self._curated_faults = None
        self._connections = None
        self._neighbour_connections = None
        self.connected_faults = None
        self._segment_dict = None
        self._multi_segment_dict = None
        self._inter_fault_connections = None
        self._trimming_gradient = trimming_gradient
        self._epsg = epsg
        self._dip_multiplier = dip_multiplier
        self._strike_multiplier = strike_multiplier
        self._sampled_dip = None
        self._trace_extensions = None
        self._additional_cuts = None
        self._excluded_cuts = None

    def add_fault(self, series: pd.Series, depth_type: str = "D90", remove_colons: bool = False,
                  tolerance: float = 100.):
        cfmFault = LeapfrogFault.from_series(series, parent_multifault=self,
                                             remove_colons=remove_colons, tolerance=tolerance,
                                             dip_choice=self.dip_choice)
        if self.smoothing_n is None:
            cfmFault.smoothed_trace = cfmFault.nztm_trace
        else:
            cfmFault.smoothed_trace = smooth_trace(cfmFault.nztm_trace, n_refinements=self.smoothing_n)
        self.faults.append(cfmFault)

    @property
    def smoothing_n(self):
        return self._smoothing_n

    @property
    def epsg(self):
        return self._epsg

    @property
    def segment_distance_tolerance(self):
        return self._segment_distance_tolerance

    @segment_distance_tolerance.setter
    def segment_distance_tolerance(self, value: float):
        self._segment_distance_tolerance = value

    @property
    def cutting_hierarchy(self):
        return self._cutting_hierarchy

    @property
    def dip_multiplier(self):
        return self._dip_multiplier

    @property
    def strike_multiplier(self):
        return self._strike_multiplier

    @cutting_hierarchy.setter
    def cutting_hierarchy(self, hierarchy: List[str]):
        for fault in self.curated_faults:
            if fault.name not in hierarchy:
                raise NameError(f"{fault.name} not included in cutting hierarchy")
        for name in hierarchy:
            if name not in self.names:
                print(f"Warning: unrecognised fault({name}) in cutting hierarchy")

        self._cutting_hierarchy = hierarchy

    def trim_cutting_hierarchy(self):
        trimmed_hierarchy = []
        for fault in self.cutting_hierarchy:
            if fault in self.names:
                trimmed_hierarchy.append(fault)
        self.cutting_hierarchy = trimmed_hierarchy

    def write_cutting_hierarchy(self, out_file: str):
        with open(out_file, "w") as out_id:
            for name in self.cutting_hierarchy:
                out_id.write(name + "\n")

    @property
    def hierarchy_dict(self):
        return self._hierarchy_dict

    @hierarchy_dict.setter
    def hierarchy_dict(self, dictionary):
        assert isinstance(dictionary, dict)
        self._hierarchy_dict = dictionary

    @property
    def curated_faults(self):
        if self._curated_faults is None:
            return self.faults
        else:
            return self._curated_faults

    @curated_faults.setter
    def curated_faults(self, fault_list: list):
        assert isinstance(fault_list, list)
        assert all([isinstance(fault, (ConnectedFaultSystem, LeapfrogFault)) for fault in fault_list])
        self._curated_faults = fault_list

    @property
    def curated_fault_dict(self):
        return {fault.name: fault for fault in self.curated_faults}

    def read_cutting_hierarchy(self, hierarchy_file: str):
        assert os.path.exists(hierarchy_file)
        with open(hierarchy_file, "r") as hfile:
            names = hfile.readlines()
            self.cutting_hierarchy = [name.strip() for name in list(names)]
        self.hierarchy_dict = {name: i for i, name in enumerate(self.cutting_hierarchy)}

    def suggest_cutting_hierarchy(self, prefix: str):
        name_ls = []
        sr_ls = []
        for fault in self.curated_faults:
            if isinstance(fault, LeapfrogFault):
                name_ls.append(fault.name)
                sr_ls.append(fault.sr_best)

            else:
                assert isinstance(fault, ConnectedFaultSystem)
                name_ls.append(fault.name)
                sr_ls.append(max([seg.sr_best for seg in fault.segments]))

        df = pd.DataFrame({"name": name_ls, "sr": sr_ls})
        df.sort_values(["sr", "name"], ascending=(False, True), inplace=True)

        out_names = list(df.name)
        out_file_name = prefix + "_suggested.csv"

        with open(out_file_name, "w") as out_id:
            for name in out_names:
                out_id.write(name + "\n")





    @property
    def names(self):
        return [fault.name for fault in self.curated_faults]
    
    @property
    def name_dict(self):
        return {fault.name: fault for fault in self.curated_faults}

    def find_connections(self, verbose: bool = True):
        """
        Find all connections between faults in the fault list using networkx
        :param verbose: print out information about individual connections
        :return:
        """
        # Total number of connections of any type
        connections = []
        # Connections between faults that are neighbours (along strike), everything else is a branch
        neighbour_connections = []
        for fault in self.faults:
            for other_fault in self.faults:
                if other_fault.name != fault.name:
                    if fault.nztm_trace.distance(other_fault.nztm_trace) <= self.segment_distance_tolerance:
                        if verbose:
                            print(f"Connection: {fault.name} and {other_fault.name}")

                        connections.append([fault.name, other_fault.name])
                        # Check if faults are neighbours
                        conditions = []
                        for p1, p2 in product([fault.end1, fault.end2], [other_fault.end1, other_fault.end2]):
                            conditions.append(p1.distance(p2) <= self.segment_distance_tolerance)
                        if any(conditions):
                            neighbour_connections.append([fault.name, other_fault.name])


        self._connections = connections
        self._neighbour_connections = neighbour_connections

        connections_set = set(tuple(sorted(conn)) for conn in connections)
        neighbour_connections_set = set(tuple(sorted(conn)) for conn in neighbour_connections)
        print(f"Found {len(connections_set)} connections")
        print(f"Found {len(neighbour_connections_set)} connections between segment ends")

    @property
    def connections(self):
        if self._connections is None:
            self.find_connections(verbose=False)
        return self._connections

    @property
    def neighbour_connections(self):
        if self._neighbour_connections is None:
            self.find_connections(verbose=False)
        return self._neighbour_connections

    def suggest_fault_systems(self, out_prefix: str):
        connected = connected_nodes(self.neighbour_connections)
        out_name = out_prefix + "_suggested.csv"
        with open(out_name, "w") as out_id:
            for connected_set in connected:
                suggested_name = suggest_combined_name(connected_set)
                out_list = [suggested_name] + sorted(list(connected_set))
                out_str = ",".join(out_list) + "\n"
                out_id.write(out_str)

    def read_fault_systems(self, fault_system_csv: str, trimming_gradient: float = 1.):
        self.connected_faults = []
        with open(fault_system_csv, "r") as in_id:
            con_data = in_id.readlines()
            for line in con_data:
                elements = line.strip().split(",")
                name = elements[0].strip()
                name = name.replace(":", "")
                segs = [element.strip() for element in elements[1:]]
                if any([seg in self.names for seg in segs]):
                    print(f"Creating connected fault system: {name} with segments: {[seg for seg in segs if seg]}")
                    cfault = ConnectedFaultSystem(overall_name=name, cfm_faults=self, segment_names=segs,
                                                  tolerance=self.segment_distance_tolerance,
                                                  trimming_gradient=trimming_gradient,
                                                  smooth_trace_refinements=self.smoothing_n)
                    self.connected_faults.append(cfault)

    def generate_curated_faults(self):
        curated_faults = []
        segments_already_included = []
        if self.connected_faults is not None:
            if len(self.connected_faults):
                curated_faults += self.connected_faults
                for fault in curated_faults:
                    segments_already_included += [seg.name for seg in fault.segments]

        else:
            print("Warning: no multi-segment faults: is this what you wanted?")

        for fault in self.faults:
            if fault.name not in segments_already_included:
                curated_faults.append(fault)
                segments_already_included.append(fault.name)

        self.curated_faults = curated_faults

    @classmethod
    def from_nz_cfm_shp(cls, filename: str, exclude_region_polygons: List[Polygon] = None, depth_type: str = "D90",
                        exclude_region_min_sr: float = 1.8, include_names: list = None, exclude_aus: bool = True,
                        exclude_zero: bool = True, sort_sr: bool = False, remove_colons: bool = False,
                        smoothing_n: Union[int, None] = 5, dip_choice: str = "pref", epsg: int = None,
                        trimming_gradient: float = 1.0, dip_multiplier: float = 1., strike_multiplier: float = 0.5,
                        check_optional_fields: bool = True):

        trimmed_fault_gdf = cls.gdf_from_nz_cfm_shp(filename=filename, exclude_region_polygons=exclude_region_polygons,
                                                    depth_type=depth_type, exclude_region_min_sr=exclude_region_min_sr,
                                                    include_names=include_names, exclude_aus=exclude_aus,
                                                    exclude_zero=exclude_zero)
        multi_fault = cls(trimmed_fault_gdf, sort_sr=sort_sr,
                          remove_colons=remove_colons, smoothing_n=smoothing_n,
                          dip_choice=dip_choice, epsg=epsg, trimming_gradient=trimming_gradient,
                          dip_multiplier=dip_multiplier, strike_multiplier=strike_multiplier,
                          check_optional_fields=check_optional_fields)

        return multi_fault

    @classmethod
    def from_shp(cls, filename: str, remove_colons: bool = False, sort_sr: bool = False, epsg: int = None,
                 dip_choice: str = "pref", trimming_gradient: float = 1.,
                 dip_multiplier: float = 1., strike_multiplier: float = 0.5, check_optional_fields: bool = False):
        assert os.path.exists(filename)
        trimmed_fault_gdf = cls.gdf_from_nz_cfm_shp(filename=filename, exclude_region_polygons=None, depth_type=None,
                                                    exclude_region_min_sr=1.8, include_names=None,
                                                    exclude_aus=False,
                                                    exclude_zero=False, )
        multi_fault = cls(trimmed_fault_gdf, sort_sr=sort_sr, remove_colons=remove_colons, epsg=epsg,
                          dip_choice=dip_choice, trimming_gradient=trimming_gradient,
                          dip_multiplier=dip_multiplier, strike_multiplier=strike_multiplier,
                          check_optional_fields=check_optional_fields)
        return multi_fault

    @property
    def segment_dict(self):
        if self._segment_dict is None:
            self.generate_segment_dicts()
        return self._segment_dict

    @property
    def multi_segment_dict(self):
        if self._multi_segment_dict is None:
            self.generate_segment_dicts()
        return self._multi_segment_dict

    def generate_segment_dicts(self):
        segment_dict = {}
        multi_segment_dict = {}
        for fault in self.faults:
            segment_dict[fault.name] = fault
        self._segment_dict = segment_dict

        for fault in self.connected_faults:
            for seg in fault.segments:
                multi_segment_dict[seg.name] = fault
        self._multi_segment_dict = multi_segment_dict

    def find_inter_fault_connections(self):
        connection_set = set()
        for con in self.connections:
            f1 = self.segment_dict[con[0]]
            f2 = self.segment_dict[con[1]]

            f1_name = f1.name if not f1.is_segment else self.multi_segment_dict[f1.name].name
            f2_name = f2.name if not f2.is_segment else self.multi_segment_dict[f2.name].name
            if f1_name != f2_name:
                connection_set.add(tuple(sorted([f1_name, f2_name], key=lambda x: x.lower())))

        self._inter_fault_connections = connection_set

    @property
    def inter_fault_connections(self):
        if self._inter_fault_connections is None:
            self.find_inter_fault_connections()
        return self._inter_fault_connections

    def find_terminations(self, fault_name: str):
        fault_index = self.cutting_hierarchy.index(fault_name)
        terminations = []
        cons = [con for con in self.inter_fault_connections if fault_name in con]
        for con in cons:
            other_fault = [a for a in list(con) if a != fault_name][0]
            other_index = self.cutting_hierarchy.index(other_fault)
            if fault_index > other_index:
                terminations.append(con)

        return terminations

    def adjust_traces_for_terminations(self, fit_distance: float = 5.e3, extend_distance: float = 20.e3, resolution: float = 1.e3):
        for fault in self.curated_faults:
            fault.adjust_trace(extend_distance=extend_distance, fit_distance=fit_distance, resolution=resolution)

    @staticmethod
    def _bearing_to_compass(bearing: float) -> str:
        """Convert a bearing in degrees (0-360) to an 8-point compass label."""
        bearing = bearing % 360
        compass_points = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        index = int((bearing + 22.5) % 360 / 45)
        return compass_points[index]

    @staticmethod
    def _build_extended_trace(trace: LineString, nearest_end: Point, other_end: Point,
                              fit_distance: float, extend_distance: float,
                              resolution: float = 1.e3) -> LineString:
        """Return what `trace` would look like after extending the `nearest_end` end.

        Mirrors the geometry produced by GenericFault.extend_trace but does not
        mutate any underlying object — used for visualising suggested extensions.
        """
        coords = np.asarray(trace.coords)
        dim = coords.shape[1]
        end_arr = np.array(nearest_end.coords)[0][:dim]
        other_end_arr = np.array(other_end.coords)[0][:dim]
        line_dists = np.linalg.norm(coords[:, :2] - end_arr[:2], axis=1)
        coords_to_fit = coords[line_dists <= fit_distance]
        local_strike = calculate_strike_direction(coords_to_fit[:, 0], coords_to_fit[:, 1])
        strike_vec = np.array([np.sin(np.radians(local_strike)),
                               np.cos(np.radians(local_strike)), 0.])[:dim]
        extend_dists = np.arange(0, extend_distance + resolution, resolution)
        if np.dot(other_end_arr - end_arr, strike_vec) < 0:
            extra = end_arr + extend_dists[:, np.newaxis] * strike_vec
        else:
            extra = end_arr - extend_dists[:, np.newaxis] * strike_vec
        if np.allclose(end_arr[:2], coords[0, :2]):
            new_coords = np.vstack([extra[::-1], coords])
        else:
            new_coords = np.vstack([coords, extra])
        return LineString(new_coords)

    def suggest_trace_extensions(self, out_file: str, fit_distance: float = 5.e3,
                                 extend_distance: float = 40.e3, proximity_threshold: float = 1.e3,
                                 geojson_out_file: str = None):
        """Suggest trace extensions for faults that terminate against higher-priority faults.

        Generates a CSV file with columns:
            fault_name, end, extend_distance, fit_distance

        The 'end' column is a compass direction (N, NE, E, SE, S, SW, W, NW)
        indicating which end of the fault trace to extend. Where both ends
        terminate, two rows are written.

        The output file can be edited by the user and then read back with
        read_trace_extensions().

        If `geojson_out_file` is provided, also writes a GeoJSON containing the
        suggested extended trace geometries (one feature per row of the CSV) so
        the suggestions can be inspected in GIS software before being applied.
        """
        rows = []
        geom_rows = []
        seen_geom_keys = set()
        for fault in self.curated_faults:
            terms = list(set(chain(*fault.find_terminations())))
            if not terms:
                continue
            cutting_faults = [self.curated_fault_dict[name] for name in terms if name != fault.name]

            if isinstance(fault, ConnectedFaultSystem):
                trace = fault.nztm_trace
            else:
                trace = fault.nztm_trace

            coords = np.array(trace.coords)
            end1 = Point(coords[0])
            end2 = Point(coords[-1])

            for nearest_end, other_end in [(end1, end2), (end2, end1)]:
                for cutting_fault in cutting_faults:
                    if nearest_end.distance(cutting_fault.nztm_trace) < proximity_threshold:
                        end_arr = np.array(nearest_end.coords)[0]
                        # Bearing from centre of fault toward the end being extended
                        centre = coords.mean(axis=0)
                        dx = end_arr[0] - centre[0]
                        dy = end_arr[1] - centre[1]
                        bearing = normalize_bearing(np.degrees(np.arctan2(dx, dy)))
                        compass = self._bearing_to_compass(bearing)
                        rows.append({
                            "fault_name": fault.name,
                            "end": compass,
                            "extend_distance": extend_distance,
                            "fit_distance": fit_distance,
                        })
                        if geojson_out_file is not None and (fault.name, compass) not in seen_geom_keys:
                            try:
                                extended = self._build_extended_trace(
                                    trace, nearest_end, other_end,
                                    fit_distance=fit_distance,
                                    extend_distance=extend_distance,
                                )
                                geom_rows.append({
                                    "fault_name": fault.name,
                                    "end": compass,
                                    "extend_distance": extend_distance,
                                    "fit_distance": fit_distance,
                                    "geometry": extended,
                                })
                                seen_geom_keys.add((fault.name, compass))
                            except Exception as e:
                                print(f"Could not build extended trace for {fault.name} ({compass}): {e}")

        df = pd.DataFrame(rows, columns=["fault_name", "end", "extend_distance", "fit_distance"])
        df.drop_duplicates(subset=["fault_name", "end"], inplace=True)
        df.to_csv(out_file, index=False)
        print(f"Wrote {len(df)} suggested trace extensions to {out_file}")

        if geojson_out_file is not None:
            if geom_rows:
                crs = f"EPSG:{self.epsg}" if self.epsg is not None else None
                geom_df = gpd.GeoDataFrame(geom_rows, geometry="geometry", crs=crs)
                geom_df.to_file(geojson_out_file, driver="GeoJSON")
                print(f"Wrote {len(geom_df)} suggested extended trace geometries to {geojson_out_file}")
            else:
                print("No extended trace geometries to write")

    def read_trace_extensions(self, csv_file: str):
        """Read a trace-extensions CSV file.

        Expected columns: fault_name, end, extend_distance, fit_distance.
        Fault names are validated against current curated fault names.
        """
        assert os.path.exists(csv_file), f"File not found: {csv_file}"
        df = pd.read_csv(csv_file)
        for col in ("fault_name", "end"):
            assert col in df.columns, f"Missing required column: {col}"
        if "extend_distance" not in df.columns:
            df["extend_distance"] = 40.e3
        if "fit_distance" not in df.columns:
            df["fit_distance"] = 5.e3

        for name in df["fault_name"]:
            if name not in self.names:
                print(f"Warning: unrecognised fault '{name}' in trace extensions file")

        self._trace_extensions = df
        print(f"Read {len(df)} trace extensions from {csv_file}")

    def read_removed_trace_extensions(self, csv_file: str):
        """Drop trace extensions listed in a removals CSV.

        Expected columns: fault_name, end. Written by review_trace_extensions.py
        as a separate list rather than by editing the curated extensions file,
        so the hand-made edited_trace_extensions.csv stays exactly as its author
        left it. Call between read_trace_extensions() and
        apply_trace_extensions().
        """
        assert os.path.exists(csv_file), f"File not found: {csv_file}"
        assert self._trace_extensions is not None, \
            "No trace extensions loaded. Call read_trace_extensions() first."
        df = pd.read_csv(csv_file)
        for col in ("fault_name", "end"):
            assert col in df.columns, f"Missing required column: {col}"

        removals = {(str(name).strip(), str(end).strip().upper())
                    for name, end in zip(df["fault_name"], df["end"])}
        keys = list(zip(self._trace_extensions["fault_name"].astype(str).str.strip(),
                        self._trace_extensions["end"].astype(str).str.strip().str.upper()))
        keep = [key not in removals for key in keys]
        dropped = len(keep) - sum(keep)
        self._trace_extensions = self._trace_extensions[keep]

        for key in removals:
            if key not in keys:
                print(f"Warning: removal '{key[0]} ({key[1]})' matches no loaded trace extension")
        print(f"Removed {dropped} trace extensions listed in {csv_file}")

    def apply_trace_extensions(self, resolution: float = 1.e3):
        """Apply trace extensions previously loaded with read_trace_extensions().

        For each row the relevant fault end is identified by compass direction
        and the trace is extended accordingly.
        """
        assert self._trace_extensions is not None, "No trace extensions loaded. Call read_trace_extensions() first."

        for _, row in self._trace_extensions.iterrows():
            fault_name = row["fault_name"]
            if fault_name not in self.names:
                continue
            fault = self.curated_fault_dict[fault_name]
            target_compass = row["end"].strip().upper()
            ext_distance = float(row["extend_distance"])
            fit_dist = float(row["fit_distance"])

            if isinstance(fault, ConnectedFaultSystem):
                trace = fault.nztm_trace
            else:
                trace = fault.nztm_trace

            coords = np.array(trace.coords)
            centre = coords.mean(axis=0)
            end1 = Point(coords[0])
            end2 = Point(coords[-1])

            # Determine which end matches the requested compass direction
            for nearest_end, other_end in [(end1, end2), (end2, end1)]:
                end_arr = np.array(nearest_end.coords)[0]
                dx = end_arr[0] - centre[0]
                dy = end_arr[1] - centre[1]
                bearing = normalize_bearing(np.degrees(np.arctan2(dx, dy)))
                compass = self._bearing_to_compass(bearing)
                if compass == target_compass:
                    if isinstance(fault, ConnectedFaultSystem):
                        nearest_seg = min(fault.segments, key=lambda s: s.nztm_trace.distance(nearest_end))
                        nearest_seg.extend_trace(nearest_end, other_end, fit_distance=fit_dist,
                                                 extend_distance=ext_distance, resolution=resolution)
                    else:
                        fault.extend_trace(nearest_end, other_end, fit_distance=fit_dist,
                                           extend_distance=ext_distance, resolution=resolution)
                    print(f"Extended {fault_name} trace at {target_compass} end by {ext_distance/1.e3:.0f} km")
                    break

    # --- Checking cut meshes against the original traces ---

    @staticmethod
    def _mesh_outline(mesh_file: str):
        """The outline (open edges) of a fault mesh.

        Returns the mesh vertices and an array of vertex-index pairs, one per
        boundary edge, or (None, None) if the mesh has no boundary.
        """
        pv_mesh = pv.read(mesh_file)
        edges = pv_mesh.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                              manifold_edges=False, non_manifold_edges=False)
        if edges.n_cells == 0:
            return None, None
        # extract_feature_edges returns two-point line cells: [2, i, j, 2, i, j, ...]
        return edges.points, edges.lines.reshape(-1, 3)[:, 1:]

    @classmethod
    def _mesh_surface_trace(cls, mesh_file: str, top_depth: float = 0., z_tolerance: float = 1.):
        """Extract the surface trace of a fault mesh as 2-D line geometry.

        The surface trace is the part of the mesh outline whose edges have both
        vertices at `top_depth` (z = 0 for these models), merged into as few
        LineStrings as possible. Returns None if the outline never reaches that
        depth, i.e. the fault's whole surface expression has been cut away.
        """
        points, line_pairs = cls._mesh_outline(mesh_file)
        if points is None:
            return None
        at_top = np.abs(points[:, 2] - top_depth) <= z_tolerance
        top_lines = line_pairs[at_top[line_pairs[:, 0]] & at_top[line_pairs[:, 1]]]
        if not len(top_lines):
            return None
        merged = linemerge([LineString([points[i, :2], points[j, :2]]) for i, j in top_lines])
        return None if merged.is_empty else merged

    @staticmethod
    def _order_edges_into_chains(line_pairs: np.ndarray):
        """Order a set of edges into chains of consecutive vertex indices.

        Consecutive edges within a chain are genuine neighbours in the outline.
        A chain runs until it closes on itself, reaches a loose end, or meets a
        vertex where more than two edges join, since past such a junction there
        is no single "next" edge.
        """
        neighbours = {}
        for a, b in line_pairs:
            neighbours.setdefault(a, []).append(b)
            neighbours.setdefault(b, []).append(a)
        walked = set()

        def walk(start, second):
            chain = [start, second]
            walked.add(frozenset((start, second)))
            while True:
                current = chain[-1]
                onward = [n for n in neighbours[current]
                          if frozenset((current, n)) not in walked]
                if len(onward) != 1:
                    break
                walked.add(frozenset((current, onward[0])))
                chain.append(onward[0])
                if onward[0] == chain[0]:
                    break
            return chain

        chains = []
        # Loose ends and junctions first, so open chains are walked in one piece
        # rather than being picked up from the middle; whatever is left over
        # after that is a closed loop and can be started anywhere.
        for vertex in sorted(neighbours, key=lambda v: len(neighbours[v]) == 2):
            for neighbour in neighbours[vertex]:
                if frozenset((vertex, neighbour)) not in walked:
                    chains.append(walk(vertex, neighbour))
        return chains

    @staticmethod
    def _kinked_segments(chain_points: np.ndarray, angle_threshold: float):
        """Find the kinks in one ordered chain of outline points.

        A segment counts as kinked if its 3-D direction differs from that of
        *both* the segments either side of it by at least `angle_threshold`
        degrees. The first and last segment of a chain have only one neighbour
        each, so they can never qualify.

        :return: (chain points with any repeats dropped, indices of the kinked
            segments within that chain, turn angle between each pair of
            consecutive segments)
        """
        # Repeated points would give a zero-length segment with no direction.
        not_repeated = np.concatenate([[True],
                                       np.linalg.norm(np.diff(chain_points, axis=0), axis=1) > 0.])
        chain_points = chain_points[not_repeated]
        if len(chain_points) < 4:
            return chain_points, np.array([], dtype=int), np.array([])

        vectors = np.diff(chain_points, axis=0)
        units = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
        turns = np.degrees(np.arccos(np.clip(np.sum(units[:-1] * units[1:], axis=1), -1., 1.)))
        # Segment i turns by turns[i - 1] from the segment before it and by
        # turns[i] from the one after it.
        kinked = np.where((turns[:-1] >= angle_threshold) & (turns[1:] >= angle_threshold))[0] + 1
        return chain_points, kinked, turns

    @staticmethod
    def _outline_plane_normal(points: np.ndarray):
        """Best-fit plane normal for a set of outline points, used to give the
        direction of a turn a consistent sign."""
        centred = points - points.mean(axis=0)
        return np.linalg.svd(centred, full_matrices=False)[2][-1]

    @staticmethod
    def _resample_chain(chain_points: np.ndarray, spacing: float):
        """Resample an ordered chain of outline points at a uniform arc length.

        Comparing chords of a fixed length is what makes a step test independent
        of how finely the outline happens to be segmented. A 1 km tread built
        from six 170 m segments turns by 0 degrees within itself, so a
        per-segment test sees nothing there and the step at its end differs from
        only one of its neighbours; measured between 1 km chords the same step
        is unmissable.
        """
        vectors = np.diff(chain_points, axis=0)
        lengths = np.linalg.norm(vectors, axis=1)
        keep = lengths > 0.
        vectors, lengths, starts = vectors[keep], lengths[keep], chain_points[:-1][keep]
        if not len(lengths):
            return None
        distance = np.concatenate([[0.], np.cumsum(lengths)])
        if distance[-1] < 3. * spacing:
            return None
        targets = np.arange(0., distance[-1], spacing)
        segment = np.clip(np.searchsorted(distance, targets, side="right") - 1, 0, len(lengths) - 1)
        fraction = (targets - distance[segment]) / lengths[segment]
        return starts[segment] + fraction[:, np.newaxis] * vectors[segment]

    @classmethod
    def _direction_reversals(cls, chain_points: np.ndarray, normal: np.ndarray,
                             spacing: float, angle_threshold: float):
        """Find where an outline turns one way and then straight back again.

        Requiring the two turns to have *opposite* sign is what separates a
        staircase from the rest of an outline: the corner where a side edge
        meets the bottom turns once, a smooth bend keeps turning the same way,
        and only a step goes back on itself. Without that test every mesh looks
        rough, because every outline has corners.

        :return: (stepped chords as LineStrings, their step amplitudes in m)
        """
        points = cls._resample_chain(chain_points, spacing)
        if points is None or len(points) < 4:
            return [], []
        vectors = np.diff(points, axis=0)
        units = vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
        turns = np.degrees(np.arccos(np.clip(np.sum(units[:-1] * units[1:], axis=1), -1., 1.)))
        signs = np.sign(np.sum(np.cross(units[:-1], units[1:]) * normal, axis=1))
        big = turns >= angle_threshold
        reversed_here = big[:-1] & big[1:] & (signs[:-1] * signs[1:] < 0.)

        segments, amplitudes = [], []
        for i in np.where(reversed_here)[0]:
            # turns[i] happens at points[i + 1] and turns[i + 1] at points[i + 2],
            # so the chord between them is the tread of the step.
            segments.append(LineString([points[i + 1], points[i + 2]]))
            amplitudes.append(spacing * np.sin(np.radians(min(turns[i], turns[i + 1]))))
        return segments, amplitudes

    @staticmethod
    def _as_2d_linestring(trace: Union[LineString, MultiLineString]):
        """Flatten a trace to a single 2-D LineString (longest part if multipart)."""
        if isinstance(trace, MultiLineString):
            trace = linemerge(trace)
        if isinstance(trace, MultiLineString):
            trace = max(trace.geoms, key=lambda g: g.length)
        return LineString([coord[:2] for coord in trace.coords])

    def _end_compass_labels(self, original_trace: LineString):
        """Label the two ends of `original_trace` the way the trace-extensions
        CSV does: by the bearing from the centre of the trace to each end.

        The labels are taken from the original (un-extended) trace because that
        is the trace suggest_trace_extensions() saw, so the labels here line up
        with the `end` column of the CSV. Doing it on the extended trace would
        not: extending an end drags the centre towards it, which can shift a
        short fault's labels by a compass point.
        """
        coords = np.array(original_trace.coords)
        centre = coords.mean(axis=0)
        out = {}
        for key, end in (("start", coords[0]), ("end", coords[-1])):
            bearing = normalize_bearing(np.degrees(np.arctan2(end[0] - centre[0], end[1] - centre[1])))
            out[key] = self._bearing_to_compass(bearing)
        return out

    def check_surface_trace_extensions(self, mesh_dir: str, suffix: str = "_cut.obj",
                                       flag_distance: float = 2.e3, lateral_tolerance: float = 500.,
                                       top_depth: float = 0., z_tolerance: float = 1.,
                                       out_file: str = None, geojson_out_file: str = None,
                                       verbose: bool = True):
        """Compare the surface trace of each cut mesh against the original fault trace.

        Trace extensions exist to grow a fault's *subsurface* area so that it
        meets its neighbours and can be cut cleanly; they are not meant to
        survive at the surface. Where a cut mesh still reaches the surface
        beyond the end of the original (un-extended) trace, either the cut did
        not happen (fix with additional_cuts) or the extension was unnecessary
        (fix by removing it from the trace-extensions CSV).

        Each cut mesh in `mesh_dir` is read back from disk, so this can be run
        on its own without re-cutting. For every fault the parts of the mesh's
        surface trace lying further than `lateral_tolerance` from the original
        trace are reported, split by which end of the trace they run off; parts
        reaching more than `flag_distance` beyond it are flagged.

        :param mesh_dir: directory holding the cut meshes
        :param suffix: appended to the fault name to make each mesh filename
        :param flag_distance: overhang (m) beyond which an end is flagged
        :param lateral_tolerance: distance (m) from the original trace within
            which surface trace is considered coincident with it; allows for
            the smoothing and spline fitting done when meshing
        :param top_depth: z value (m) of the top of the meshes
        :param z_tolerance: tolerance (m) on `top_depth`
        :param out_file: optional CSV to write the results to
        :param geojson_out_file: optional GeoJSON of the extra trace geometry
        :param verbose: print a summary of the flagged faults
        :return: DataFrame of results, worst overhang first
        """
        assert os.path.isdir(mesh_dir), f"Directory not found: {mesh_dir}"
        extensions = self._trace_extensions
        rows = []
        missing_meshes = []
        no_surface_trace = []
        n_checked = 0

        for fault in self.curated_faults:
            mesh_file = os.path.join(mesh_dir, f"{fault.name}{suffix}")
            if not os.path.exists(mesh_file):
                missing_meshes.append(fault.name)
                continue
            if fault.original_nztm_trace is None:
                print(f"Warning: no original trace stored for {fault.name}; skipping")
                continue

            n_checked += 1
            surface_trace = self._mesh_surface_trace(mesh_file, top_depth=top_depth,
                                                     z_tolerance=z_tolerance)
            if surface_trace is None:
                no_surface_trace.append(fault.name)
                continue

            original = self._as_2d_linestring(fault.original_nztm_trace)
            compass = self._end_compass_labels(original)

            # Everything further from the original trace than lateral_tolerance is "extra".
            extra = surface_trace.difference(original.buffer(lateral_tolerance))
            if extra.is_empty:
                continue
            parts = list(extra.geoms) if hasattr(extra, "geoms") else [extra]

            # Group the extra pieces by the end of the original trace they run off.
            by_end = {}
            for part in parts:
                if not isinstance(part, LineString) or part.length == 0.:
                    continue
                coords = np.array(part.coords)
                distances = [original.distance(Point(xy)) for xy in coords]
                furthest = int(np.argmax(distances))
                position = original.project(Point(coords[furthest]))
                if position <= 0.:
                    end = compass["start"]
                elif position >= original.length:
                    end = compass["end"]
                else:
                    # Not off either end: an along-strike bulge or a stray fragment.
                    end = "side"
                record = by_end.setdefault(end, {"max": 0., "length": 0., "parts": []})
                record["max"] = max(record["max"], max(distances))
                record["length"] += part.length
                record["parts"].append(part)

            for end, record in by_end.items():
                if extensions is not None and end != "side":
                    matching = extensions[(extensions["fault_name"] == fault.name) &
                                          (extensions["end"].str.strip().str.upper() == end)]
                else:
                    matching = None
                rows.append({
                    "fault_name": fault.name,
                    "end": end,
                    "max_extra_m": record["max"],
                    "extra_length_m": record["length"],
                    "flagged": record["max"] > flag_distance,
                    "has_extension": bool(matching is not None and len(matching)),
                    "extend_distance": (float(matching["extend_distance"].iloc[0])
                                        if matching is not None and len(matching) else np.nan),
                    "geometry": linemerge(record["parts"]) if len(record["parts"]) > 1 else record["parts"][0],
                })

        columns = ["fault_name", "end", "max_extra_m", "extra_length_m", "flagged",
                   "has_extension", "extend_distance", "geometry"]
        df = pd.DataFrame(rows, columns=columns)
        df.sort_values("max_extra_m", ascending=False, inplace=True, ignore_index=True)

        # These are inspection-only outputs, so a file left open in QGIS should
        # warn rather than throw away the results of the whole check.
        if out_file is not None:
            try:
                df.drop(columns="geometry").to_csv(out_file, index=False)
                print(f"Wrote {len(df)} surface-trace checks to {out_file}")
            except PermissionError as e:
                print(f"Could not write {out_file} (open in another program?): {e}")
        if geojson_out_file is not None:
            if len(df):
                crs = f"EPSG:{self.epsg}" if self.epsg is not None else None
                try:
                    gpd.GeoDataFrame(df, geometry="geometry", crs=crs).to_file(geojson_out_file,
                                                                              driver="GeoJSON")
                    print(f"Wrote {len(df)} extra surface trace geometries to {geojson_out_file}")
                except PermissionError as e:
                    print(f"Could not write {geojson_out_file} (open in QGIS?): {e}")
            else:
                print("No extra surface trace geometries to write")

        if verbose:
            flagged = df[df["flagged"]]
            print(f"\nChecked surface traces of {n_checked} cut meshes against the original traces "
                  f"({len(df)} ends with >{lateral_tolerance:.0f} m of extra trace, "
                  f"{len(flagged)} over the {flag_distance/1.e3:.1f} km threshold)")
            if missing_meshes:
                print(f"  {len(missing_meshes)} curated faults have no mesh in {mesh_dir}")
            if no_surface_trace:
                print(f"  {len(no_surface_trace)} meshes reach the surface nowhere: "
                      f"{', '.join(no_surface_trace[:5])}"
                      f"{' ...' if len(no_surface_trace) > 5 else ''}")
            for _, row in flagged.iterrows():
                fix = "remove/shorten the trace extension, or add a cut" if row["has_extension"] \
                    else "check the cutting hierarchy"
                print(f"  {row['fault_name']} ({row['end']}): {row['max_extra_m']/1.e3:.1f} km beyond "
                      f"the original trace, {row['extra_length_m']/1.e3:.1f} km of extra trace -- {fix}")

        return df

    def check_jagged_cuts(self, mesh_dir: str, suffix: str = "_cut.obj",
                          angle_threshold: float = 20., min_kinked_segments: int = 3,
                          step_scales: tuple = (500., 1000.), step_angle: float = 30.,
                          min_reversals: int = 3,
                          top_depth: float = 0., z_tolerance: float = 1.,
                          out_file: str = None, geojson_out_file: str = None,
                          verbose: bool = True):
        """Flag cut meshes whose subsurface outline is jagged.

        A clean cut leaves a smooth outline. A cut that has gone wrong leaves a
        rough one, and roughness comes in two forms that need measuring
        differently:

        *Kinks* are single segments that differ in 3-D direction from **both**
        their neighbours by at least `angle_threshold`. This catches isolated
        spikes and the sliver segments a bad clip leaves behind. A fault with
        `min_kinked_segments` or more is flagged.

        *Steps* are the staircases a cut leaves when it follows whole triangle
        edges instead of the true intersection. These are invisible to the kink
        test, because each tread and riser is several collinear segments: the
        turn inside a tread is zero, and the segment at a corner differs from
        only one of its neighbours. They are found instead by resampling the
        outline at each of `step_scales` and counting *direction reversals* --
        consecutive turns of at least `step_angle` that go opposite ways. The
        scale with the most reversals is reported, and `min_reversals` or more
        flags the fault.

        Only the subsurface part of the outline is examined. The surface trace
        is dropped first, because a trace faithfully follows the mapped fault
        and is expected to bend sharply; use check_surface_trace_extensions()
        for that part of the outline.

        :param mesh_dir: directory holding the cut meshes
        :param suffix: appended to the fault name to make each mesh filename
        :param angle_threshold: direction change (degrees) counting as a kink
        :param min_kinked_segments: number of kinks needed to flag a fault
        :param step_scales: chord lengths (m) at which to look for steps
        :param step_angle: direction change (degrees) counting as a step turn
        :param min_reversals: number of reversals needed to flag a fault
        :param top_depth: z value (m) of the top of the meshes
        :param z_tolerance: tolerance (m) on `top_depth`
        :param out_file: optional CSV to write the results to
        :param geojson_out_file: optional GeoJSON of the kinked and stepped parts
        :param verbose: print a summary of the flagged faults
        :return: DataFrame of every mesh with a kink or a step, worst first
        """
        assert os.path.isdir(mesh_dir), f"Directory not found: {mesh_dir}"
        rows = []
        n_checked = 0

        for fault in self.curated_faults:
            mesh_file = os.path.join(mesh_dir, f"{fault.name}{suffix}")
            if not os.path.exists(mesh_file):
                continue
            n_checked += 1
            points, line_pairs = self._mesh_outline(mesh_file)
            if points is None:
                continue

            # Drop the surface trace, leaving the subsurface outline. Chains
            # then break where the trace used to be, so the segments dropping
            # away from the surface are chain ends and are never counted as
            # kinks against the trace they hang off.
            at_top = np.abs(points[:, 2] - top_depth) <= z_tolerance
            subsurface = line_pairs[~(at_top[line_pairs[:, 0]] & at_top[line_pairs[:, 1]])]
            if not len(subsurface):
                continue

            chains = [points[np.array(chain)]
                      for chain in self._order_edges_into_chains(subsurface)]
            normal = self._outline_plane_normal(points)

            segments = []
            turn_angles = []
            longest_run = 0
            for chain_points in chains:
                chain_points, kinked, turns = self._kinked_segments(chain_points, angle_threshold)
                if not len(kinked):
                    continue
                for i in kinked:
                    segments.append(LineString([chain_points[i], chain_points[i + 1]]))
                    turn_angles.append(max(turns[i - 1], turns[i]))
                # Kinks that run together are a serrated stretch of outline
                # rather than isolated spikes, so record the longest such run.
                run = np.split(kinked, np.where(np.diff(kinked) != 1)[0] + 1)
                longest_run = max(longest_run, max(len(part) for part in run))

            # Steps, at whichever of the scales shows them most clearly.
            step_segments, step_amplitudes, step_scale = [], [], np.nan
            for scale in step_scales:
                scale_segments, scale_amplitudes = [], []
                for chain_points in chains:
                    found, amplitudes = self._direction_reversals(chain_points, normal, scale,
                                                                  step_angle)
                    scale_segments.extend(found)
                    scale_amplitudes.extend(amplitudes)
                if len(scale_segments) > len(step_segments):
                    step_segments, step_amplitudes, step_scale = (scale_segments,
                                                                  scale_amplitudes, scale)

            geometry = segments + step_segments
            if not geometry:
                continue
            rows.append({
                "fault_name": fault.name,
                "n_kinked_segments": len(segments),
                "longest_run": longest_run,
                "max_turn_deg": max(turn_angles) if turn_angles else 0.,
                "n_reversals": len(step_segments),
                "reversal_scale_m": step_scale,
                "step_amplitude_m": (float(np.median(step_amplitudes)) if step_amplitudes
                                     else np.nan),
                "flagged_kinks": len(segments) >= min_kinked_segments,
                "flagged_steps": len(step_segments) >= min_reversals,
                "flagged": (len(segments) >= min_kinked_segments
                            or len(step_segments) >= min_reversals),
                "geometry": MultiLineString(geometry) if len(geometry) > 1 else geometry[0],
            })

        columns = ["fault_name", "n_kinked_segments", "longest_run", "max_turn_deg",
                   "n_reversals", "reversal_scale_m", "step_amplitude_m",
                   "flagged_kinks", "flagged_steps", "flagged", "geometry"]
        df = pd.DataFrame(rows, columns=columns)
        # Steps first: they are the defect that usually means the cut went wrong.
        df.sort_values(["n_reversals", "n_kinked_segments"], ascending=False, inplace=True,
                       ignore_index=True)

        # Inspection-only outputs, so a file left open elsewhere should warn
        # rather than throw away the results of the whole check.
        if out_file is not None:
            try:
                df.drop(columns="geometry").to_csv(out_file, index=False)
                print(f"Wrote {len(df)} jagged-cut checks to {out_file}")
            except PermissionError as e:
                print(f"Could not write {out_file} (open in another program?): {e}")
        if geojson_out_file is not None:
            if len(df):
                crs = f"EPSG:{self.epsg}" if self.epsg is not None else None
                try:
                    gpd.GeoDataFrame(df, geometry="geometry", crs=crs).to_file(geojson_out_file,
                                                                              driver="GeoJSON")
                    print(f"Wrote kinked outline segments for {len(df)} faults to {geojson_out_file}")
                except PermissionError as e:
                    print(f"Could not write {geojson_out_file} (open in QGIS?): {e}")
            else:
                print("No kinked outline segments to write")

        if verbose:
            flagged = df[df["flagged"]]
            print(f"\nChecked the subsurface outline of {n_checked} cut meshes for jaggedness "
                  f"({len(flagged)} flagged of {len(df)} with any roughness: "
                  f"{int(df['flagged_steps'].sum())} stepped, "
                  f"{int(df['flagged_kinks'].sum())} kinked)")
            for _, row in flagged.iterrows():
                detail = []
                if row["flagged_steps"]:
                    detail.append(f"{row['n_reversals']} reversals at "
                                  f"{row['reversal_scale_m']:.0f} m, steps of "
                                  f"{row['step_amplitude_m']:.0f} m")
                if row["flagged_kinks"]:
                    detail.append(f"{row['n_kinked_segments']} kinked segments "
                                  f"(longest run {row['longest_run']}, worst turn "
                                  f"{row['max_turn_deg']:.0f} deg)")
                print(f"  {row['fault_name']}: {'; '.join(detail)}")

        return df

    def candidate_cutters(self, fault_name: str):
        """Every fault that could have cut `fault_name`: those ranked above it
        in the cutting hierarchy, plus any forced onto it by additional_cuts
        (which apply even where the named cutter ranks lower, or is absent
        from the hierarchy altogether)."""
        if fault_name in self.cutting_hierarchy:
            higher = self.cutting_hierarchy[:self.cutting_hierarchy.index(fault_name)]
        else:
            higher = []
        forced = [cut for (cut_fault, cut) in self.additional_cuts
                  if cut_fault == fault_name and cut not in higher]
        return higher + forced

    def attribute_jagged_cuts(self, jagged: pd.DataFrame, uncut_mesh_dir: str,
                              suffix: str = "_depth_contours.obj", max_distance: float = 1.5e3,
                              depth_surface: pv.PolyData = None, verify_cuts: bool = True,
                              actual_cuts=None,
                              cut_threshold: float = 0.7, min_cut_distance: float = 10.e3,
                              bottom_depth: float = -30000., out_file: str = None,
                              geojson_out_file: str = None, verbose: bool = True):
        """Work out which cut left each jagged fault's outline serrated.

        Takes the output of check_jagged_cuts() and, for every kinked segment,
        finds the nearest fault surface that was entitled to cut that fault
        (see candidate_cutters()). A kink sitting on a cutter's surface was
        left by that cut. Kinks near neither, but close to the base depth
        surface, are attributed to it; the rest come back "unattributed" and
        are worth a look, since they point at something other than a cut --
        awkward contours, or a surface self-intersecting.

        The uncut Stage 1 surfaces are used for the cutters, because those are
        the meshes the cutting loop starts from and they span the whole
        intersection even where the cutter was later trimmed itself.

        Proximity alone will happily blame a fault that merely passes nearby,
        so each attribution is checked and reported in the `cut_confirmed`
        column. Pass `actual_cuts` -- the record the cutting loop kept of what
        it really cut (see FaultMesh.cut_against) -- and the column is exact.
        Without it, `verify_cuts` falls back to re-running the model's own cut
        decision; that test is run without the higher-priority context that can
        veto a cut, which only makes it more willing to say yes, so a False
        still means the real run cannot have made that cut either.

        Attributions that come back False are worth reading the other way round
        -- the outline is following a surface that never cut it, which may be a
        cut that should have been made and was not.

        :param jagged: DataFrame from check_jagged_cuts(), typically the
            flagged rows only
        :param uncut_mesh_dir: directory of the uncut (Stage 1) surfaces
        :param suffix: appended to the fault name to make each mesh filename
        :param max_distance: how close (m) a kink must lie to a surface to be
            attributed to it
        :param depth_surface: the base depth surface, if kinks should also be
            attributed to the depth trim
        :param verify_cuts: check each attribution when `actual_cuts` is not given
        :param actual_cuts: the cuts really made, as {fault_name: [cutter, ...]}
            or the path to a two-column cuts_applied.csv; makes `cut_confirmed`
            exact and skips the re-test entirely
        :param cut_threshold: `threshold` passed to decide_whether_to_cut
        :param min_cut_distance: `min_distance` passed to decide_whether_to_cut
        :param bottom_depth: `bottom_depth` passed to decide_whether_to_cut
        :param out_file: optional CSV to write the results to
        :param geojson_out_file: optional GeoJSON of the attributed segments
        :param verbose: print the attribution per fault
        :return: DataFrame with a row per (jagged fault, cutter)
        """
        assert os.path.isdir(uncut_mesh_dir), f"Directory not found: {uncut_mesh_dir}"

        # Read every candidate surface once; the cutters are drawn from the
        # whole hierarchy, so there is little to be saved by being lazy here.
        vertices = {}
        for fault in self.curated_faults:
            path = os.path.join(uncut_mesh_dir, f"{fault.name}{suffix}")
            if os.path.exists(path):
                vertices[fault.name] = meshio.read(path).points
        if verbose:
            print(f"Read {len(vertices)} uncut surfaces from {uncut_mesh_dir} for attribution")
        extents = {name: np.array([points[:, 0].min(), points[:, 0].max(),
                                   points[:, 1].min(), points[:, 1].max()])
                   for name, points in vertices.items()}

        # KD-trees are only built for surfaces that survive the extent test.
        trees = {}

        def tree_for(name):
            if name not in trees:
                trees[name] = cKDTree(vertices[name])
            return trees[name]

        depth_tree = None
        if depth_surface is not None:
            depth_points = np.asarray(depth_surface.points)
            depth_tree = cKDTree(depth_points[:, :2])

        rows = []
        for _, jagged_row in jagged.iterrows():
            fault_name = jagged_row["fault_name"]
            geometry = jagged_row["geometry"]
            segments = list(geometry.geoms) if hasattr(geometry, "geoms") else [geometry]
            midpoints = np.array([np.asarray(seg.coords).mean(axis=0) for seg in segments])

            kink_extent = np.array([midpoints[:, 0].min(), midpoints[:, 0].max(),
                                    midpoints[:, 1].min(), midpoints[:, 1].max()])
            nearest_cutter = np.full(len(segments), None, dtype=object)
            nearest_distance = np.full(len(segments), np.inf)
            for cutter in self.candidate_cutters(fault_name):
                if cutter not in vertices:
                    continue
                xmin, xmax, ymin, ymax = extents[cutter]
                if (xmin > kink_extent[1] + max_distance or xmax < kink_extent[0] - max_distance or
                        ymin > kink_extent[3] + max_distance or ymax < kink_extent[2] - max_distance):
                    continue
                distances, _ = tree_for(cutter).query(midpoints)
                closer = distances < nearest_distance
                nearest_distance[closer] = distances[closer]
                nearest_cutter[closer] = cutter

            if depth_tree is not None:
                _, index = depth_tree.query(midpoints[:, :2])
                depth_distance = np.abs(midpoints[:, 2] - depth_points[index, 2])
            else:
                depth_distance = np.full(len(segments), np.inf)

            by_cutter = {}
            for i, segment in enumerate(segments):
                if nearest_distance[i] <= max_distance:
                    cutter, distance = nearest_cutter[i], nearest_distance[i]
                elif depth_distance[i] <= max_distance:
                    cutter, distance = "base depth surface", depth_distance[i]
                else:
                    cutter, distance = "unattributed", nearest_distance[i]
                record = by_cutter.setdefault(cutter, {"segments": [], "distances": []})
                record["segments"].append(segment)
                record["distances"].append(distance)

            for cutter, record in sorted(by_cutter.items(),
                                         key=lambda item: -len(item[1]["segments"])):
                rows.append({
                    "fault_name": fault_name,
                    "cutter": cutter,
                    "n_segments": len(record["segments"]),
                    "median_distance_m": float(np.median(record["distances"])),
                    "geometry": (MultiLineString(record["segments"])
                                 if len(record["segments"]) > 1 else record["segments"][0]),
                })

        if actual_cuts is not None:
            self._confirm_from_record(rows, actual_cuts)
        elif verify_cuts:
            self._verify_attributions(rows, uncut_mesh_dir, suffix, cut_threshold,
                                      min_cut_distance, bottom_depth)

        columns = ["fault_name", "cutter", "n_segments", "median_distance_m"]
        if verify_cuts or actual_cuts is not None:
            columns.append("cut_confirmed")
        columns.append("geometry")
        df = pd.DataFrame(rows, columns=columns)

        # Inspection-only outputs, so a file left open elsewhere should warn
        # rather than throw away the results of the whole run.
        if out_file is not None:
            try:
                df.drop(columns="geometry").to_csv(out_file, index=False)
                print(f"Wrote {len(df)} jagged-cut attributions to {out_file}")
            except PermissionError as e:
                print(f"Could not write {out_file} (open in another program?): {e}")
        if geojson_out_file is not None:
            if len(df):
                crs = f"EPSG:{self.epsg}" if self.epsg is not None else None
                try:
                    gpd.GeoDataFrame(df, geometry="geometry", crs=crs).to_file(geojson_out_file,
                                                                              driver="GeoJSON")
                    print(f"Wrote {len(df)} attributed kink groups to {geojson_out_file}")
                except PermissionError as e:
                    print(f"Could not write {geojson_out_file} (open in QGIS?): {e}")
            else:
                print("No attributed kinks to write")

        if verbose and len(df):
            unattributed = df[df["cutter"] == "unattributed"]["n_segments"].sum()
            print(f"\nAttributed the kinks on {df['fault_name'].nunique()} faults "
                  f"({df['n_segments'].sum() - unattributed} of {df['n_segments'].sum()} "
                  f"segments matched to a cut)")
            if verify_cuts or actual_cuts is not None:
                # Explicit identity test: the column also holds pd.NA, which
                # makes an == comparison ambiguous.
                unconfirmed = df["cut_confirmed"].apply(lambda value: value is False).sum()
                print(f"  {unconfirmed} attributions are to a fault that did not cut this "
                      f"one (cut_confirmed False) -- a surface being followed without a cut")
            for fault_name, group in df.groupby("fault_name", sort=False):
                blame = ", ".join(
                    f"{row['cutter']} x{row['n_segments']}"
                    f"{'?' if verify_cuts and row['cut_confirmed'] is False else ''}"
                    for _, row in group.iterrows())
                print(f"  {fault_name}: {blame}")

        return df

    @staticmethod
    def _confirm_from_record(rows: list, actual_cuts):
        """Fill in `cut_confirmed` from the record of cuts really made.

        `actual_cuts` is either {fault_name: [cutter, ...]} or the path to a
        two-column cuts_applied.csv written by the cutting loop. Nothing is
        re-tested here: a cut either happened or it did not.
        """
        if not isinstance(actual_cuts, dict):
            frame = pd.read_csv(actual_cuts)
            actual_cuts = {name: list(group["cutter"])
                           for name, group in frame.groupby("fault_name")}
        for row in rows:
            if row["cutter"] in ("unattributed", "base depth surface"):
                row["cut_confirmed"] = pd.NA
            else:
                row["cut_confirmed"] = row["cutter"] in actual_cuts.get(row["fault_name"], [])

    def _verify_attributions(self, rows: list, uncut_mesh_dir: str, suffix: str,
                             cut_threshold: float, min_cut_distance: float,
                             bottom_depth: float):
        """Put each attributed (fault, cutter) pair back through the model's own
        cut decision, adding a `cut_confirmed` entry to every row in place."""
        meshes = {}

        def mesh_for(name):
            if name not in meshes:
                path = os.path.join(uncut_mesh_dir, f"{name}{suffix}")
                meshes[name] = FaultMesh.from_file(path) if os.path.exists(path) else None
            return meshes[name]

        decisions = {}
        for row in rows:
            fault_name, cutter = row["fault_name"], row["cutter"]
            if cutter in ("unattributed", "base depth surface"):
                row["cut_confirmed"] = pd.NA
                continue
            if (fault_name, cutter) not in decisions:
                decision = self.should_cut(fault_name, cutter)
                if decision is None:
                    target_mesh, cutting_mesh = mesh_for(fault_name), mesh_for(cutter)
                    if target_mesh is None or cutting_mesh is None:
                        decision = pd.NA
                    else:
                        try:
                            # higher_meshes is left out on purpose: it can only
                            # veto a cut, so omitting it keeps a False answer
                            # meaning "the real run cannot have cut this either".
                            decision = bool(target_mesh.decide_whether_to_cut(
                                cutting_mesh, threshold=cut_threshold,
                                min_distance=min_cut_distance, higher_meshes=None,
                                bottom_depth=bottom_depth, fancy_cutting=True))
                        except Exception as e:
                            print(f"  Could not test cut of {fault_name} by {cutter}: {e}")
                            decision = pd.NA
                decisions[(fault_name, cutter)] = decision
            row["cut_confirmed"] = decisions[(fault_name, cutter)]

    # --- Additional and excluded cuts ---

    @property
    def additional_cuts(self):
        """Set of (fault_name, fault_name) tuples for cuts that should always be made."""
        if self._additional_cuts is None:
            return set()
        return self._additional_cuts

    @property
    def excluded_cuts(self):
        """Set of (fault_name, fault_name) tuples for cuts that should never be made."""
        if self._excluded_cuts is None:
            return set()
        return self._excluded_cuts

    def read_additional_cuts(self, csv_file: str):
        """Read a two-column CSV listing pairs of faults that should always be cut.

        Columns: fault_to_cut, cutting_fault
        Fault names are validated against current curated fault names.
        """
        assert os.path.exists(csv_file), f"File not found: {csv_file}"
        df = pd.read_csv(csv_file, header=None, names=["fault_to_cut", "cutting_fault"])
        cuts = set()
        for _, row in df.iterrows():
            name1 = row["fault_to_cut"].strip()
            name2 = row["cutting_fault"].strip()
            for name in (name1, name2):
                if name not in self.names:
                    print(f"Warning: unrecognised fault '{name}' in additional cuts file")
            cuts.add((name1, name2))

        self._additional_cuts = cuts
        print(f"Read {len(cuts)} additional cuts from {csv_file}")

    def read_excluded_cuts(self, csv_file: str):
        """Read a two-column CSV listing pairs of faults that should NOT be cut.

        Columns: fault_to_cut, cutting_fault
        Fault names are validated against current curated fault names.
        """
        assert os.path.exists(csv_file), f"File not found: {csv_file}"
        df = pd.read_csv(csv_file, header=None, names=["fault_to_cut", "cutting_fault"])
        cuts = set()
        for _, row in df.iterrows():
            name1 = row["fault_to_cut"].strip()
            name2 = row["cutting_fault"].strip()
            for name in (name1, name2):
                if name not in self.names:
                    print(f"Warning: unrecognised fault '{name}' in excluded cuts file")
            cuts.add((name1, name2))

        self._excluded_cuts = cuts
        print(f"Read {len(cuts)} excluded cuts from {csv_file}")

    def should_cut(self, fault_name: str, cutting_fault_name: str) -> bool:
        """Check additional/excluded cut lists and return True/False/None.

        Returns True if the pair is in additional_cuts, False if in excluded_cuts,
        or None if neither list applies (fall through to decide_whether_to_cut).
        """
        if (fault_name, cutting_fault_name) in self.excluded_cuts:
            print(f"Skipping cut of {fault_name} by {cutting_fault_name} (excluded cuts list).")
            return False
        if (fault_name, cutting_fault_name) in self.additional_cuts:
            print(f"Forcing cut of {fault_name} by {cutting_fault_name} (additional cuts list).")
            return True
        return None

    def to_opensha_xml(self, exclude_subduction: bool = True, buffer_width: float = 5000.,
                       write_buffers: bool = True, subduction_names: tuple = ("hikurangi", "puysegur")):
        """
        Write out XML in OpenSHA format
        :param exclude_subduction: Do not include subduction zones from CFM
        :return:
        """
        assert self.faults
        assert isinstance(exclude_subduction, bool)
        # Base XML element
        opensha_element = ElemTree.Element("OpenSHA")
        # Fault model sub element
        fm_element = ElemTree.Element("FaultModel")
        opensha_element.append(fm_element)

        i = 0
        for fault in self.faults:
            # Identify subduction zone sources to include (only if exclude_subduction is True)
            exclude_condition = all([exclude_subduction, any([name in fault.name.lower()
                                                              for name in subduction_names])])
            # Add XML for fault
            if not exclude_condition:
                fm_element.append(fault.to_xml(section_id=i, buffer_width=buffer_width, write_buffers=write_buffers))
                i += 1

        # Awkward way of getting the xml file to be written in a way that's easy to read.
        # elmstr = ElemTree.tostring(opensha_element, encoding="UTF-8", xml_declaration=True)
        elmstr = ElemTree.tostring(opensha_element, encoding="UTF-8")
        xml_dom = minidom.parseString(elmstr)
        pretty_xml_str = xml_dom.toprettyxml(indent="  ", encoding="utf-8")

        return pretty_xml_str
    
    def to_quads_mesh(self, sampled_dip: bool = True, depth_multiplier: float = 0.8, file_name: str = None):
        """
        Generate a mesh of quads from the fault traces.
        :param sampled_dip: If True, use the sampled dip value for the fault.
        :param depth_multiplier: Multiplier for the depth to control the spacing of the mesh.
        :return: Meshio mesh object.
        """
        assert self.faults, "No faults to generate mesh from"
        
        meshes = []
        for fault in self.faults:
            vertices, quads = fault.to_quads(sampled=sampled_dip, depth_multiplier=depth_multiplier)
            mesh = meshio.Mesh(points=vertices, cells={"quad": quads})
            
            meshes.append(mesh)

        merged_meshes = pv.merge([pv.from_meshio(m) for m in meshes])
        if file_name is not None:
            if not file_name.endswith(".vtk"):
                file_name += ".vtk"
            merged_meshes.save(file_name)

        return merged_meshes

    def combine_meshes_to_vtk(self, file_name: str = None, only_faults: List[str] = None,
                              include_connected_segments: bool = True):
        """
        Merge each curated fault's surface mesh into a single VTK file.

        Each fault contributes its `.mesh.mesh` (a meshio Mesh on the
        attached FaultMesh). A `fault_id` cell array distinguishes faults
        in the merged output; the id → name mapping is stored as field
        data under `fault_names`. For `ConnectedFaultSystem` faults whose
        top-level mesh isn't set, falls back to merging individual segment
        meshes when `include_connected_segments` is True.

        Args:
            file_name (str): Optional output path. Defaults to no save.
                `.vtk` extension is appended if missing.
            only_faults (List[str]): Optional restricted list of fault
                names to include.
            include_connected_segments (bool): If True, when a
                ConnectedFaultSystem has no combined mesh, merge in any
                meshes attached to its individual segments.

        Returns:
            pv.PolyData (or UnstructuredGrid): The merged mesh.
        """
        pieces = []
        fault_names: List[str] = []
        for fault in self.curated_faults:
            if only_faults is not None and fault.name not in only_faults:
                continue

            fault_meshes = []
            top_mesh = getattr(fault, "_mesh", None)
            if top_mesh is not None and getattr(top_mesh, "mesh", None) is not None:
                fault_meshes.append(top_mesh)
            elif include_connected_segments and isinstance(fault, ConnectedFaultSystem):
                for seg in fault.segments:
                    seg_mesh = getattr(seg, "_mesh", None)
                    if seg_mesh is not None and getattr(seg_mesh, "mesh", None) is not None:
                        fault_meshes.append(seg_mesh)

            if not fault_meshes:
                continue

            fault_id = len(fault_names)
            for fm in fault_meshes:
                piece = pv.from_meshio(fm.mesh)
                piece.cell_data["fault_id"] = np.full(piece.n_cells, fault_id, dtype=np.int32)
                pieces.append(piece)
            fault_names.append(fault.name)

        if not pieces:
            raise RuntimeError("No fault meshes available to combine.")

        merged = pv.merge(pieces)
        merged.field_data["fault_names"] = np.array(fault_names)

        if file_name is not None:
            if not file_name.endswith(".vtk"):
                file_name += ".vtk"
            merged.save(file_name)

        return merged


class LeapfrogFault(GenericFault):
    """
    Represents either a whole fault (for simple faults) or one segment. Behaviours is slightly
    """
    def __init__(self, parent_multifault: LeapfrogMultiFault = None, smoothing: int = 5,
                 trimming_gradient: float = 1.0, segment_distance_tolerance: float = 100.):
        self._end1 = None
        self._end2 = None
        self._neighbouring_segments = None
        self._neighbour_dict = {}
        self._neighbour_angle_dict = {}

        super(LeapfrogFault, self).__init__(None)
        self._parent = parent_multifault

        self.connections = []
        self._is_segment = False
        self.parent_connected_fault = None
        self._smoothing = smoothing
        self._trimming_gradient = trimming_gradient
        self._segment_distance_tolerance = segment_distance_tolerance
        self._contours = None
        self._smoothed_trace = None
        self._footprint = None
        self._sampled_dip = None
        self._mesh = None

    @property
    def is_segment(self):
        """
        Records whether instance is a segment of a larger multi-segment fault like the Alpine Fault.
        :return:
        """
        return self._is_segment

    @property
    def smoothing(self):
        """
        n value to use in Chaikin's corner-cutting algorithm.
        :return:
        """
        return self._smoothing

    @property
    def dip_multiplier(self):
        return self.parent._dip_multiplier

    @property
    def strike_multiplier(self):
        return self.parent._strike_multiplier

    @property
    def trimming_gradient(self):
        """
        Factor that controls how much the ends of segment contours of a multi-segment fault are
        shortened to allow
        :return:
        """
        return self.parent._trimming_gradient

    @property
    def neighbouring_segments(self):
        return self._neighbouring_segments

    @neighbouring_segments.setter
    def neighbouring_segments(self, segment_list: list):
        assert isinstance(segment_list, list)
        if len(segment_list) > 2:
            raise ValueError(f"Too many ({len(segment_list):d}) neighbours supplied for segment {self.name:s}.\n"
                             f"Only two neighbours are allowed. Please turn extra neighbour(s) into separate faults.\n"
                             f"Neighbours supplied: {[neighbour.name for neighbour in segment_list]}")
        if len(segment_list) == 0:
            self._is_segment = False
        else:
            self._is_segment = True

        self._neighbouring_segments = segment_list

        for seg in segment_list:
            if seg.nztm_trace.distance(self.end1) <= self.segment_distance_tolerance:
                self._neighbour_dict[self.end1.coords[0]] = seg
                self._neighbour_angle_dict[self.end1.coords[0]] = self.neighbour_angle(seg)
            elif seg.nztm_trace.distance(self.end2) <= self.segment_distance_tolerance:
                self._neighbour_dict[self.end2.coords[0]] = seg
                self._neighbour_angle_dict[self.end2.coords[0]] = self.neighbour_angle(seg)
            else:
                raise ValueError(f"Orphan neighbour segment for {self.name}: {seg.name}")

    @property
    def neighbour_angle_dict(self):
        return self._neighbour_angle_dict

    @property
    def neighbour_dict(self):
        return self._neighbour_dict

    @property
    def strike(self):
        return normalize_bearing(self.dip_dir - 90.)

    @property
    def along_strike_vector(self):
        strike_rad = np.radians(self.strike)
        return np.array([np.sin(strike_rad), np.cos(strike_rad), 0.])

    @property
    def across_strike_vector(self):
        dip_dir_rad = np.radians(self.dip_dir)
        return np.array([np.sin(dip_dir_rad), np.cos(dip_dir_rad), 0.])

    @property
    def smoothed_trace(self):
        return self._smoothed_trace

    @smoothed_trace.setter
    def smoothed_trace(self, trace: LineString):
        assert isinstance(trace, LineString)
        self._smoothed_trace = trace

    @property
    def parent(self):
        """
        Return LeapfrogMultiFault instance that this fault is part of.
        :return:
        """
        return self._parent
    
    @property
    def sampled_dip(self):
        """
        Dip value sampled from the fault trace. Used for calculating dip in depth contours.
        :return:
        """
        if self._sampled_dip is None:
            self.generate_sampled_dip(from_parent=False)
        return self._sampled_dip
    
    @sampled_dip.setter
    def sampled_dip(self, dip: float):
        """
        Set the sampled dip value.
        :param dip: Dip value in degrees.
        """
        assert isinstance(dip, (int, float))
        assert 0. <= dip <= 180., "Dip must be between 0 and 180 degrees"
        self._sampled_dip = dip

    def generate_sampled_dip(self, from_parent = True):
        """
        Generate a sampled dip value from the min and max dip
        :param from_parent: If True, use the parent multi-fault's sampled dip value.
        :return:
        """
        if from_parent and self.parent is not None:
            self.sampled_dip = self.dip_min + (self.dip_max - self.dip_min) * self.parent.sampled_dip
        else:
            self.sampled_dip = self.dip_min + (self.dip_max - self.dip_min) * np.random.rand()

    @property
    def down_dip_vector_sampled(self):
        """
        Calculated from dip and dip direction
        """
        assert self.sampled_dip is not None
        if self.dip_dir is None:
            # Assume vertical
            return np.array([0., 0., -1])
        else:
            z = np.sin(np.radians(self.sampled_dip))
            x, y = np.cos(np.radians(self.sampled_dip)) * np.array([np.sin(np.radians(self.dip_dir)),
                                                                np.cos(np.radians(self.dip_dir))])
        return np.array([x, y, -z])

    def depth_contour(self, depth: float, smoothing: bool = True, km=False, sr_and_rake: bool = False, dummy_segment_length: float = 100.):
        """
        Generate contour of fault surface at depth below surface
        :param depth: In metres, upwards is positive
        :param smoothing: N for use with Chaikin's corner cutting
        :param km: If True, divide depth by 1000
        :return: LineString or MultiLineString representing contour
        """
        if depth <= 0:
            shift = depth / self.down_dip_vector[-1]
        else:
            shift = (-1 * depth) / self.down_dip_vector[-1]
        if km:
            shift *= 1000.

        xo, yo, zo = shift * self.down_dip_vector

        if smoothing:
            assert self.smoothing is not None
            if self.is_segment:
                smoothed_contour = self.smoothed_trace

            else:
                self.smoothed_trace = smooth_trace(self.nztm_trace, self.smoothing)
                smoothed_contour = self.smoothed_trace


        else:
            smoothed_contour = self.nztm_trace

        shifted_contour = translate(smoothed_contour, xo, yo, zo)

        if self.is_segment:

            contour_e1 = Point(shifted_contour.coords[0])
            contour_e2 = Point(shifted_contour.coords[-1])

            flipping_conditions = [all([self.end1.x > self.end2.x, contour_e1.x < contour_e2.x]),
                                   all([self.end1.x < self.end2.x, contour_e1.x > contour_e2.x]),
                                   all([self.end1.x == self.end2.x, self.end1.y > self.end2.y,
                                        contour_e1.y < contour_e2.y]),
                                   all([self.end1.x == self.end2.x, self.end1.y < self.end2.y,
                                        contour_e1.y > contour_e2.y])]

            if any(flipping_conditions):
                contour_e2 = Point(shifted_contour.coords[0])
                contour_e1 = Point(shifted_contour.coords[-1])

            if len(self.neighbouring_segments) > 1:
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                e1_box_intersection = e1_box.boundary.intersection(shifted_contour)
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                e2_box_intersection = e2_box.boundary.intersection(shifted_contour)

                if (not e1_box.intersects(e2_box)
                        and isinstance(e1_box_intersection, Point)
                        and isinstance(e2_box_intersection, Point)):
                    trimmed_contour = cut_line_between_two_points(shifted_contour, [e1_box_intersection,
                                                                                     e2_box_intersection])
                else:
                    trimmed_contour = None

            elif self.end1.coords[0] in self.neighbour_dict.keys():
                e1_box = self.end_clipping_box(self.end1, depth, gradient_adjustment=self.trimming_gradient)
                if all([point.within(e1_box) for point in (contour_e1, contour_e2)]):
                    trimmed_contour = None
                else:
                    e1_box_intersection = e1_box.boundary.intersection(shifted_contour)

                    if isinstance(e1_box_intersection, Point):
                        split_line = cut_line_at_point(shifted_contour, e1_box_intersection)
                        if split_line[0].distance(contour_e1) == 0.:
                            trimmed_contour = split_line[1]
                        elif split_line[1].distance(contour_e1) == 0.:
                            trimmed_contour = split_line[0]
                        else:
                            print("neither intersects")
                            trimmed_contour = None
                    else:
                        
                        trimmed_contour = None
                if trimmed_contour is None:
                    contour_e1_array = np.array(contour_e1.coords[0])
                    contour_e2_array = np.array(contour_e2.coords[0])
                    e1_e2_vector = contour_e2_array - contour_e1_array
                    if np.dot(e1_e2_vector, self.along_strike_vector) >= 0:
                        trimmed_contour = LineString([contour_e2, contour_e2_array + self.along_strike_vector * dummy_segment_length])
                    else:
                        trimmed_contour = LineString([contour_e2, contour_e2_array - self.along_strike_vector * dummy_segment_length])

                    

            else:  # self.end2 in self.neighbour_dict.keys()
                e2_box = self.end_clipping_box(self.end2, depth, gradient_adjustment=self.trimming_gradient)
                if all([point.within(e2_box) for point in (contour_e1, contour_e2)]):
                    trimmed_contour = None
                else:
                    e2_box_intersection = e2_box.boundary.intersection(shifted_contour)
                    if isinstance(e2_box_intersection, Point):
                        split_line = cut_line_at_point(shifted_contour, e2_box_intersection)
                        if split_line[0].distance(contour_e2) == 0.:
                            trimmed_contour = split_line[1]
                        elif split_line[1].distance(contour_e2) == 0.:
                            trimmed_contour = split_line[0]
                        else:
                            print("neither intersects")
                            trimmed_contour = None
                    else:
                        trimmed_contour = None
                
                if trimmed_contour is None:
                    contour_e1_array = np.array(contour_e1.coords[0])
                    contour_e2_array = np.array(contour_e2.coords[0])
                    e1_e2_vector = contour_e2_array - contour_e1_array
                    if np.dot(e1_e2_vector, self.along_strike_vector) >= 0:
                        trimmed_contour = LineString([contour_e1, contour_e1_array - self.along_strike_vector * dummy_segment_length])
                    else:
                        trimmed_contour = LineString([contour_e1, contour_e1_array + self.along_strike_vector * dummy_segment_length])


            # if trimmed_contour is not None:
            #     if abs(trimmed_contour.length - smoothed_contour.length) < distance_tolerance:
            #         trimmed_contour = None

        else:
            trimmed_contour = shifted_contour
        if sr_and_rake:
            sr = self.sr_best
            rake = self.rake_best
            if trimmed_contour is not None:
                return trimmed_contour, sr, rake
            else:
                return None, None, None
        else:
            return trimmed_contour
        

    def generate_depth_contours(self, depths: Union[np.ndarray, List[float]], smoothing: bool = False,
                       km: bool = False, trace: bool = False):
        if trace:
            contours = [self.nztm_trace] + [self.depth_contour(depth, smoothing, km) for depth in depths]
        else:
            contours = [self.depth_contour(depth, smoothing, km) for depth in depths]

        if max(depths) > 0:
            depths *= -1

        if self.parent.epsg is not None:
            self.contours =  gpd.GeoDataFrame({"depth": depths}, geometry=contours,crs=self.parent.epsg)
        else:
            self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours)
    
    def depth_contour_mesh(self, depths: Union[np.ndarray, List[float]], smoothing: bool = False,
                               km: bool = False, trace: bool = False, point_spacing: float = 1000., output_spacing: float = 1000.,
                               mesh_name: Union[str, None] = None, mesh_format: str = "vtk", check_mesh: bool = False,
                               check_strike_dip: bool = True):
        if trace:
            contours = [self.nztm_trace] + [self.depth_contour(depth, smoothing, km) for depth in depths]
        else:
            contours = [self.depth_contour(depth, smoothing, km) for depth in depths]

        if max(depths) > 0:
            depths *= -1

        if self.parent.epsg is not None:
            self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours, crs=self.parent.epsg)
        else:
            self.contours = gpd.GeoDataFrame({"depth": depths}, geometry=contours)

        contour_spline = spline_fit_contours(self.contours, point_spacing=point_spacing, output_spacing=output_spacing)

        if mesh_name is None:
            mesh = triangulate_contours(contour_spline, mesh_format=mesh_format, check_mesh=check_mesh,
                                        check_strike_dip=check_strike_dip)
        else:
            mesh = triangulate_contours(contour_spline, mesh_name=mesh_name, mesh_format=mesh_format,
                                        check_mesh=check_mesh, check_strike_dip=check_strike_dip)


        return mesh
    
    def mesh_fault_surface(self, resolution: float = 1000., spline_resolution: float = 100., plane_fitting_eps: float = 1.0e-5, check_mesh: bool = False, check_strike_dip: bool = False):
        if self.contours is None:
            raise ValueError("No contours to mesh. Please run generate_depth_contours() first.")

        contour_spline = spline_fit_contours(self.contours, point_spacing=spline_resolution, output_spacing=resolution)

        mesh = triangulate_contours(contour_spline, mesh_format="vtk", check_mesh=check_mesh,
                                    check_strike_dip=check_strike_dip)
        return mesh
        

    @property
    def nztm_trace(self):
        return self._nztm_trace
    
    @property
    def nztm_trace_coords(self):
        """
        Return the coordinates of the nztm trace as a numpy array.
        :return: numpy array of coordinates
        """
        return np.array(self.nztm_trace.coords)

    @nztm_trace.setter
    def nztm_trace(self, trace: LineString):
        assert isinstance(trace, (LineString, MultiLineString))
        if isinstance(trace, MultiLineString):
            try:
                trace = merge_multiple_nearly_adjacent_segments(list(trace.geoms))
            except Exception as e:
                print(f"Error merging segments: {e}. Taking first segment only.")
                trace = list(trace.geoms)[0]

        if trace.has_z:
            new_trace = LineString([(xi, yi, 0.) for xi, yi, _ in trace.coords])
        else:
            new_trace = LineString([(xi, yi, 0.) for xi, yi in trace.coords])
        self._nztm_trace = new_trace
        new_coord_array = np.array(new_trace.coords)
        self._end1 = Point(new_coord_array[0])
        self._end2 = Point(new_coord_array[-1])

    @property
    def nztm_trace_geoseries(self):
        if self.parent.epsg is not None:
            return gpd.GeoSeries(self.nztm_trace, crs=self.parent.epsg)
        else:
            return gpd.GeoSeries(self.nztm_trace)
        
    @property
    def nztm_trace_array(self):
        return np.array(self.nztm_trace.coords)
    
    @property
    def original_nztm_trace(self):
        return self._original_nztm_trace
    
    @original_nztm_trace.setter
    def original_nztm_trace(self, trace: LineString):
        assert isinstance(trace, (LineString, MultiLineString))
        if isinstance(trace, MultiLineString):
            trace = list(trace.geoms)[0]

        if trace.has_z:
            new_trace = LineString([(xi, yi, 0.) for xi, yi, _ in trace.coords])
        else:
            new_trace = LineString([(xi, yi, 0.) for xi, yi in trace.coords])
        self._original_nztm_trace = new_trace

    @property
    def original_nztm_trace_array(self):
        return np.array(self.original_nztm_trace.coords)

    @property
    def end1(self):
        return self._end1

    @property
    def end2(self):
        return self._end2

    def clipping_box(self, centre_point: Point, along_half_width: float, across_half_width: float = 100000.):
        point_array = np.array(centre_point.coords)
        across_shift = across_half_width * self.across_strike_vector
        along_shift = along_half_width * self.along_strike_vector
        out_array = point_array + np.array([across_shift + along_shift,
                                            across_shift - along_shift,
                                            - 1 * (across_shift + along_shift),
                                            along_shift - across_shift])
        return Polygon(out_array)


    def end_clipping_box(self, end_i: Point, depth: float, gradient_adjustment: float = 1.,
                         across_half_width: float = 100000.):
        end_angle = self.neighbour_angle_dict[end_i.coords[0]]
        end_width = gradient_adjustment * np.tan(np.radians(end_angle)) * (depth / np.sin(np.radians(self.dip_best)))

        return self.clipping_box(end_i, end_width, across_half_width=across_half_width)

    def neighbour_angle(self, neighbour: LeapfrogFault):
        strike_diff = smallest_difference(self.dip_dir, neighbour.dip_dir)
        if strike_diff > 90.:
            strike_diff = 180. - strike_diff
            dip_diff = 180 - self.dip_best - neighbour.dip_best
            if dip_diff > 90.:
                print(f"Warning: {self.name} and {neighbour.name}: shallow dips in opposite directions."
                       "Are you sure you want to connect?")
        else:
            dip_diff = abs(neighbour.dip_best - self.dip_best)
        return max([self.dip_multiplier * dip_diff, self.strike_multiplier * strike_diff])

    @property
    def segment_distance_tolerance(self):
        if self.parent is not None:
            return self.parent.segment_distance_tolerance
        else:
            return self._segment_distance_tolerance

    @property
    def footprint(self):
        if self._footprint is None:
            self.calculate_footprint()
        return self._footprint

    @property
    def footprint_geoseries(self):
        if self.parent.epsg is not None:
            return gpd.GeoSeries(self.footprint, crs=self.parent.epsg)
        else:
            return gpd.GeoSeries(self.footprint)
    @property
    def footprint_linestring(self):
        return LineString(self.footprint.exterior.coords)

    def calculate_footprint(self, smoothed: bool = True, buffer: float = 5000.):
        if smoothed:
            trace = self.smoothed_trace
        else:
            trace = self.nztm_trace

        buffer_top_offset = self.across_strike_vector * -1 * buffer
        shifted_top = translate(trace, *buffer_top_offset)

        bottom_trace = list(self.contours.geometry)[-1]
        buffer_bottom_offset = self.across_strike_vector * buffer
        shifted_bottom = translate(bottom_trace, *buffer_bottom_offset)

        buffer_combined = MultiLineString([shifted_top, shifted_bottom])
        self._footprint = buffer_combined.minimum_rotated_rectangle

    def find_terminations(self):
        return self.parent.find_terminations(self.name)

    def adjust_footprint(self):
        terms = list(set(chain(*self.find_terminations())))
        if terms:
            cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
            for fault in cutting_faults:
                for nearest_end, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
                    if nearest_end.distance(fault.nztm_trace) < 1.e3:
                        if isinstance(fault, ConnectedFaultSystem):
                            closest_seg = min(fault.segments, key=lambda x: x.nztm_trace.distance(nearest_end))
                        else:
                            closest_seg = fault

                        self.extend_footprint(nearest_end, other_end, closest_seg)

            # fp_to_merge = [self.footprint] + [fault.footprint for fault in cutting_faults]
            # merged_footprints = unary_union(fp_to_merge)
            # cutting_lines = []
            # for end_i, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
            #     if not any([end_i.distance(fault.nztm_trace) < self.segment_distance_tolerance for fault in cutting_faults]):
            #         cutting_lines.append((LineString([np.array(end_i) + self.across_strike_vector * line_length,
            #                                           np.array(end_i) - self.across_strike_vector * line_length]),
            #                               other_end))
            # if len(cutting_lines):
            #     if len(cutting_lines) > 1:
            #         print(f"{self.name}: more than one cutting line, choosing first...")
            #     splitter, other_end = cutting_lines[0]
            #     split_footprint = split(merged_footprints, splitter)
            #     kept_polys = [poly for poly in list(split_footprint) if other_end.within(poly)]
            #     if len(kept_polys) > 1:
            #         print(f"{self.name}: more than one cut polygon, choosing first...")
            #     self._footprint = kept_polys[0]
            #
            # else:
            #     self._footprint = merged_footprints

    def extend_footprint(self, end_i: Point, other_end: Point, other_segment: LeapfrogFault,
                         deepest_contour_depth: float = 30.e3, search_line_length: float = 1.5e5,
                         buffer_size: float = 5.e3, fall_back_distance: float = 40.3):
        """

        :param end_i: End to extend
        :param other_end: Other end of segment
        :param other_segment: Other
        :param deepest_contour_depth:
        :param search_line_length:
        :param buffer_size:
        :param fall_back_distance:
        :return:
        """
        # Find strike direction
        diff_vector = np.array(other_end.coords) - np.array(end_i.coords)
        if np.dot(diff_vector, self.along_strike_vector) > 0:
            strike_direction = -1 * self.along_strike_vector
        else:
            strike_direction = self.along_strike_vector


        if self.is_segment:
            bottom_trace = None
            contour_depth = deepest_contour_depth
            while bottom_trace is None:
                bottom_trace = self.depth_contour(contour_depth)
                contour_depth -= 2000.
        else:
            bottom_trace = list(self.contours.geometry)[-1]
        bot1 = np.array(bottom_trace.coords)[0]
        bot2 = np.array(bottom_trace.coords)[-1]

        def distance_along_strike(point: np.array, reference: np.array, strike_vector: np.array):
            diff = point - np.array(reference)
            return np.dot(diff, strike_vector)

        bot_i = bot1 if (distance_along_strike(bot1, np.array(end_i.coords), strike_direction) >
                         distance_along_strike(bot2, np.array(end_i.coords), strike_direction)) else bot2

        search_line = LineString([bot_i,
                                  bot_i + search_line_length * strike_direction])

        if other_segment.is_segment:
            other_fault = other_segment.parent_connected_fault
        else:
            other_fault = other_segment

        other_contour = other_fault.depth_contour(deepest_contour_depth)
        if other_contour.intersects(bottom_trace):
            return

        else:
            trace_intersection = search_line.intersection(other_fault.nztm_trace)
            if isinstance(trace_intersection, MultiPoint):
                trace_intersection = min(list(trace_intersection.geoms), key=lambda x: x.distance(Point(bot_i)))
            contour_intersection = search_line.intersection(other_contour)
            if isinstance(contour_intersection, MultiPoint):
                contour_intersection = min(list(contour_intersection.geoms), key=lambda x: x.distance(Point(bot_i)))

            if any([a.is_empty for a in [trace_intersection, contour_intersection]]):
                if all([a.is_empty for a in [trace_intersection, contour_intersection]]):
                    corner_point = Point(bot_i + fall_back_distance * strike_direction)
                elif not trace_intersection.is_empty:
                    corner_point = trace_intersection
                else:
                    corner_point = contour_intersection
            else:
                corner_point = max([contour_intersection, trace_intersection], key=lambda x: x.distance(Point(bot_i)))

            triangle = Polygon(np.vstack([bot_i, np.array(corner_point.coords), np.array(end_i.coords)])).buffer(buffer_size,
                                                                                        cap_style=2)


            if self.is_segment:
                new_boundary = unary_union([self.parent_connected_fault.footprint, triangle])
                self.parent_connected_fault._footprint = new_boundary
            else:
                new_boundary = unary_union([self.footprint, triangle])
                self._footprint = new_boundary

            return
        
    def extend_trace(self, end_i: Point, other_end: Point,
                     fit_distance: float = 5.e3, extend_distance: float = 20.3, resolution: float = 1.e3):
        end_arr = np.array(end_i.coords)[0]
        other_end_arr = np.array(other_end.coords)[0]
        line_dists = np.linalg.norm(self.nztm_trace_coords - end_arr, axis=1)
        coords_to_fit = self.nztm_trace_coords[line_dists <= fit_distance]
        local_strike_direction = calculate_strike_direction(coords_to_fit[:, 0], coords_to_fit[:, 1])
        local_strike_vector = np.array([np.sin(np.radians(local_strike_direction)), np.cos(np.radians(local_strike_direction)), 0.])
        extend_dists = np.arange(0, extend_distance + resolution, resolution)

        if np.dot(other_end_arr - end_arr, local_strike_vector) < 0:
            # Extend the trace in the direction of the local strike vector
            extra_coords = end_arr + extend_dists[:, np.newaxis] * local_strike_vector
        else:
            # Extend the trace in the opposite direction
            extra_coords = end_arr - extend_dists[:, np.newaxis] * local_strike_vector

        original_nztm = self.nztm_trace_coords.copy()
        if np.allclose(end_arr, original_nztm[0]):
            new_nztm_coords = np.vstack([extra_coords[-1::-1], original_nztm])
        else:
            assert np.allclose(end_arr, original_nztm[-1])
            new_nztm_coords = np.vstack([original_nztm, extra_coords])

        self.nztm_trace = LineString(new_nztm_coords)

    def adjust_trace(self, extend_distance: float = 20.e3, fit_distance: float = 5.e3, resolution: float = 1.e3):
        # Adjust the trace by extending and fitting it
        terms = list(set(chain(*self.find_terminations())))
        if terms:
            cutting_faults = [self.parent.curated_fault_dict[name] for name in terms if name is not self.name]
            for fault in cutting_faults:
                for nearest_end, other_end in zip([self.end1, self.end2], [self.end2, self.end1]):
                    if nearest_end.distance(fault.nztm_trace) < 1.e3:
                        self.extend_trace(nearest_end, other_end, fit_distance=fit_distance, extend_distance=extend_distance, resolution=resolution)

    def to_xml(self, section_id: int, buffer_width: float = 5000., write_buffers: bool = True):
        # Unique fault identifier
        tag_name = "i{:d}".format(section_id)
        # Metadata
        attribute_dic = {"sectionId": "{:d}".format(section_id),
                         "sectionName": self.name,
                         "aveLongTermSlipRate": "{:.2f}".format(self.sr_best),
                         "slipRateStdDev": "{:.2f}".format(self.sr_sigma),
                         "aveDip": "{:.1f}".format(self.dip_best),
                         "aveRake": "{:.1f}".format(self.rake_to_opensha(self.rake_best)),
                         "aveUpperDepth": "0.0",
                         "aveLowerDepth": "{:.1f}".format(self.depth_best),
                         "aseismicSlipFactor": "0.0",
                         "couplingCoeff": "1.0",
                         "dipDirection": "{:.1f}".format(self.dip_dir),
                         "parentSectionId": "-1",
                         "connector": "false",
                         "domainNo": f"{self.dom_num:d}",
                         "domainName": f"{self.dom_name}"
                         }

        # Initialize XML element
        fault_element = ElemTree.Element(tag_name, attrib=attribute_dic)
        # Add sub element for fault trace
        trace_element = fault_trace_xml(self.wgs_trace, self.name)
        fault_element.append(trace_element)
        if write_buffers:
            # Add sub element for FZ buffer
            polygon_element = fault_polygon_xml(self.combined_buffer_polygon(buffer_width), self.name)
            fault_element.append(polygon_element)

    def to_quads(self, sampled: bool = True, depth_multiplier: float = 0.8) -> dict:
        """
        Extract quads and vertices for writing using meshio.
        :param sampled: If True, use the sampled dip vector.
        :param depth_multiplier: Multiplier for the depth to scale the bottom trace.
        :return: vertices and quads as numpy arrays.
        """
        top_trace = np.array(self.nztm_trace.coords)
        if sampled:
            bottom_trace = top_trace + self.down_dip_vector_sampled * depth_multiplier * self.depth_best * 1000. * -1./ self.down_dip_vector_sampled[-1]
        else:
            bottom_trace = top_trace + self.down_dip_vector * depth_multiplier * self.depth_best * 1000. * -1./ self.down_dip_vector[-1]

        vertices = np.vstack([top_trace, bottom_trace])

        quad_list = []
        for (top_i, bottom_i) in zip(range(len(top_trace)-1), range(len(bottom_trace)-1)):
            quad_list.append([top_i, top_i + 1, len(top_trace) + bottom_i + 1, len(top_trace) + bottom_i])

        quad_array = np.array(quad_list, dtype=int)

        return vertices, quad_array
    
    @property
    def mesh(self):
        return self._mesh
    
    @mesh.setter
    def mesh(self, mesh: FaultMesh):
        assert isinstance(mesh, FaultMesh)
        self._mesh = mesh


            


