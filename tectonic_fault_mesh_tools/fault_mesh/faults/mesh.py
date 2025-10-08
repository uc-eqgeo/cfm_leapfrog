import meshio
import numpy as np
import os
from fault_mesh.utilities.meshing import fit_plane_to_points, get_strike_dip_from_normal
from itertools import combinations
import pyvista as pv

class FaultMesh:
    """
    A class representing a fault mesh.

    Attributes:
        vertices (list): A list of vertices in the mesh.
        faces (list): A list of faces in the mesh.
    """

    def __init__(self):
        """
        Initializes a FaultMesh instance.

        Args:
            vertices (list, optional): A list of vertices. Defaults to an empty list.
            faces (list, optional): A list of faces. Defaults to an empty list.
        """
        self.mesh = None
        self.vertices = None
        self.triangles = None
        self.tri_dict = None
        self.vertex_dict = None
        self.plane = None
        self._plane_normal = None
        self._plane_origin = None
        self._plane_along_strike = None
        self._plane_down_dip = None
        self._plane_strike = None
        self._plane_dip = None
        self._resolved_vertex_x = None
        self._resolved_vertex_y = None
        self._resolved_vertex_z = None
        self.name = None


    def __repr__(self):
        return f"FaultMesh(name={self.name}, vertices={len(self.vertices)}, triangles={len(self.triangles)})"

    @classmethod
    def from_file(cls, file_path):
        """
        Loads a mesh from a file.

        Args:
            file_path (str): The path to the mesh file.
            
        Returns:
            FaultMesh: A FaultMesh instance with the loaded data.
        """
        assert os.path.isfile(file_path), f"File {file_path} does not exist."

        mesh = meshio.read(file_path)
        fault_mesh = cls()
        fault_mesh.mesh = mesh
        fault_mesh.vertices = mesh.points
        fault_mesh.triangles = mesh.cells_dict["triangle"]
        fault_mesh.tri_dict = {i: tri for i, tri in enumerate(fault_mesh.triangles)}
        fault_mesh.vertex_dict = {i: vertex for i, vertex in enumerate(fault_mesh.vertices)}
        fault_mesh.name = os.path.splitext(os.path.basename(file_path))[0]
        return fault_mesh

    def fit_plane(self, eps=1.0e-5):
        """
        Fits a plane to the mesh vertices and returns the plane normal and origin.

        Args:
            eps (float, optional): Tolerance for plane fitting. Defaults to 1.0e-5.

        Returns:
            tuple: A tuple containing the plane normal and origin.
        """
        assert self.vertices is not None, "Mesh vertices are not defined."

        plane_normal, plane_origin = fit_plane_to_points(self.vertices, eps=eps)

        self._plane_normal = plane_normal
        self._plane_origin = plane_origin
        self._plane_strike, self._plane_dip = get_strike_dip_from_normal(plane_normal)
        strike_rad = np.radians(self._plane_strike)
        self._plane_along_strike = np.array([np.sin(strike_rad), np.cos(strike_rad), 0])
        down_dip = np.cross(plane_normal, self._plane_along_strike)
        if down_dip[2] > 0:
            down_dip = -down_dip
        self._plane_down_dip = down_dip / np.linalg.norm(down_dip)

        return plane_normal, plane_origin
    
    @property
    def plane_normal(self):
        if self._plane_normal is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_normal

    @property
    def plane_origin(self):
        if self._plane_origin is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_origin
    
    @property
    def plane_strike(self):
        if self._plane_strike is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_strike

    @property
    def plane_dip(self):
        if self._plane_dip is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_dip
    
    @property
    def plane_along_strike(self):
        if self._plane_along_strike is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_along_strike
    
    @property
    def plane_down_dip(self):
        if self._plane_down_dip is None:
            self._plane_normal, self._plane_origin = self.fit_plane()
        return self._plane_down_dip
    
    @property
    def resolved_vertex_x(self):
        if self._resolved_vertex_x is None:
            assert self.vertices is not None, "Mesh vertices are not defined."
            self._resolved_vertex_x = np.dot(self.vertices - self.plane_origin, self.plane_along_strike)
        return self._resolved_vertex_x

    @property
    def resolved_vertex_y(self):
        if self._resolved_vertex_y is None:
            assert self.vertices is not None, "Mesh vertices are not defined."
            self._resolved_vertex_y = np.dot(self.vertices - self.plane_origin, self.plane_down_dip)
        return self._resolved_vertex_y
    
    @property
    def resolved_vertex_z(self):
        if self._resolved_vertex_z is None:
            assert self.vertices is not None, "Mesh vertices are not defined."
            self._resolved_vertex_z = np.dot(self.vertices - self.plane_origin, self.plane_normal)
        return self._resolved_vertex_z
    
    def find_all_outside_edges(self):
        triangle_edges = np.array([np.sort(np.array(list(combinations(tri, 2))), axis=1) for tri in self.triangles])
        unique_edges, edge_counts = np.unique(triangle_edges.reshape(-1, 2), axis=0, return_counts=True)
        outside_edges = unique_edges[edge_counts == 1]
        return outside_edges

    def find_all_outside_vertex_indices(self):
        outside_edges = self.find_all_outside_edges()
        outside_vertex_indices = np.unique(outside_edges.reshape(-1, 2))
        return outside_vertex_indices

    def find_all_outside_vertices(self):
        outside_vertex_indices = self.find_all_outside_vertex_indices()
        return self.vertices[outside_vertex_indices]
    
    def find_vertex_indices_constant_depth(self, bottom_value: float, tolerance=10.0):
        outside_vertex_indices = self.find_all_outside_vertex_indices()
        bottom_vertex_indices = np.where(np.abs(self.vertices[:, -1][outside_vertex_indices] - bottom_value) <= tolerance)[0]
        return outside_vertex_indices[bottom_vertex_indices]

    def find_edge_triangles(self, edge):
        """
        Finds triangles that contain the given edge.

        Args:
            edge (array-like): A 2-element array-like object representing the edge vertices.

        Returns:
            list: A list of triangle indices that contain the edge.
        """
        edge = np.sort(edge)
        triangle_indices = []
        for i, tri in self.tri_dict.items():
            if edge[0] in tri and edge[1] in tri:
                triangle_indices.append(i)
        return triangle_indices
    
    def find_left_right_vertices(self, top_depth: float = 0., bottom_depth: float = -30000., tolerance: float = 10., min_separation=5.e3):
        """Finds the left and right vertices of the mesh at a given depth range.

        :param top_depth: The top depth to consider, defaults to 0.
        :type top_depth: float, optional
        :param bottom_depth: The bottom depth to consider, defaults to -20000.
        :type bottom_depth: float, optional
        :param tolerance: _description_, defaults to 10.
        :type tolerance: float, optional
        :param min_separation: _description_, defaults to 5.e3
        :type min_separation: _type_, optional
        """

        top_vertex_indices = self.find_vertex_indices_constant_depth(top_depth, tolerance=tolerance)
        bottom_vertex_indices = self.find_vertex_indices_constant_depth(bottom_depth, tolerance=tolerance)
        all_edge_vertices = self.find_all_outside_vertex_indices()
        side_vertex_indices = np.setdiff1d(all_edge_vertices, np.concatenate([top_vertex_indices, bottom_vertex_indices]))
        side_vertex_x = self.resolved_vertex_x[side_vertex_indices]
        sorted_side_vertex_order = np.argsort(side_vertex_x)
        sorted_side_vertex_indices = side_vertex_indices[sorted_side_vertex_order]
        sorted_side_vertex_x = side_vertex_x[sorted_side_vertex_order]
        biggest_gap = np.max(np.diff(sorted_side_vertex_x))
        if biggest_gap < min_separation:
            print(f"No significant gap found in side vertices: {self.name}")
            return None, None
        gap_index = np.argmax(np.diff(sorted_side_vertex_x))
        left_indices = sorted_side_vertex_indices[:gap_index + 1]
        right_indices = sorted_side_vertex_indices[gap_index + 1:]
        return left_indices, right_indices
    
    def find_top_bottom_left_right_vertices(self, top_depth: float = 0., bottom_depth: float = -30000., tolerance: float = 10., min_separation=5.e3):
        """Finds the top, bottom, left, and right vertices of the mesh at a given depth range.

        :param top_depth: The top depth to consider, defaults to 0.
        :type top_depth: float, optional
        :param bottom_depth: The bottom depth to consider, defaults to -20000.
        :type bottom_depth: float, optional
        :param tolerance: _description_, defaults to 10.
        :type tolerance: float, optional
        :param min_separation: _description_, defaults to 5.e3
        :type min_separation: _type_, optional
        """

        top_vertex_indices = self.find_vertex_indices_constant_depth(top_depth, tolerance=tolerance)
        bottom_vertex_indices = self.find_vertex_indices_constant_depth(bottom_depth, tolerance=tolerance)
        assert len(top_vertex_indices) > 0, f"No top vertices found at depth {top_depth} +/- {tolerance}"
        if len(bottom_vertex_indices) == 0:
            print(f"No bottom vertices found at depth {bottom_depth} +/- {tolerance}, using deepest vertices instead.")
            bottom_vertex_indices = np.array([self.find_all_outside_vertex_indices()[np.argmin(self.vertices[self.find_all_outside_vertex_indices(), -1])]])
        all_edge_vertices = self.find_all_outside_vertex_indices()
        side_vertex_indices = np.setdiff1d(all_edge_vertices, np.concatenate([top_vertex_indices, bottom_vertex_indices]))
        side_vertex_x = self.resolved_vertex_x[side_vertex_indices]
        sorted_side_vertex_order = np.argsort(side_vertex_x)
        sorted_side_vertex_indices = side_vertex_indices[sorted_side_vertex_order]
        sorted_side_vertex_x = side_vertex_x[sorted_side_vertex_order]
        biggest_gap = np.max(np.diff(sorted_side_vertex_x))
        if biggest_gap < min_separation:
            print(f"No significant gap found in side vertices: {self.name}")
            return None, None, None, None
        gap_index = np.argmax(np.diff(sorted_side_vertex_x))
        left_indices = sorted_side_vertex_indices[:gap_index + 1]
        right_indices = sorted_side_vertex_indices[gap_index + 1:]
        return top_vertex_indices, bottom_vertex_indices, left_indices, right_indices   
    
    def find_top_bottom_left_right_edges(self, top_depth: float = 0., bottom_depth: float = -30000., tolerance: float = 10., min_separation=5.e3):
        """Finds the top, bottom, left, and right edges of the mesh at a given depth range.

        :param top_depth: The top depth to consider, defaults to 0.
        :type top_depth: float, optional
        :param bottom_depth: The bottom depth to consider, defaults to -20000.
        :type bottom_depth: float, optional
        :param tolerance: _description_, defaults to 10.
        :type tolerance: float, optional
        :param min_separation: _description_, defaults to 5.e3
        :type min_separation: _type_, optional
        """

        top_vertex_indices, bottom_vertex_indices, left_vertex_indices, right_vertex_indices = self.find_top_bottom_left_right_vertices(top_depth=top_depth, bottom_depth=bottom_depth, tolerance=tolerance, min_separation=min_separation)
        if top_vertex_indices is None or bottom_vertex_indices is None or left_vertex_indices is None or right_vertex_indices is None:
            return None, None, None, None
        all_edges = self.find_all_outside_edges()
        top_edges = []
        bottom_edges = []
        left_edges = []
        right_edges = []

        for edge in all_edges:
            if edge[0] in top_vertex_indices and edge[1] in top_vertex_indices:
                top_edges.append(edge)
            elif edge[0] in bottom_vertex_indices and edge[1] in bottom_vertex_indices:
                bottom_edges.append(edge)
            if any([vert_i in left_vertex_indices for vert_i in edge]):
                    left_edges.append(edge)
            elif any([vert_i in right_vertex_indices for vert_i in edge]):
                    right_edges.append(edge)

        return np.array(top_edges), np.array(bottom_edges), np.array(left_edges), np.array(right_edges)
    
    def find_top_bottom_left_right_triangles(self, top_depth: float = 0., bottom_depth: float = -30000., tolerance: float = 10., min_separation=5.e3):
        """Finds the top, bottom, left, and right triangles of the mesh at a given depth range.

        :param top_depth: The top depth to consider, defaults to 0.
        :type top_depth: float, optional
        :param bottom_depth: The bottom depth to consider, defaults to -20000.
        :type bottom_depth: float, optional
        :param tolerance: _description_, defaults to 10.
        :type tolerance: float, optional
        :param min_separation: _description_, defaults to 5.e3
        :type min_separation: _type_, optional
        """

        top_edges, bottom_edges, left_edges, right_edges = self.find_top_bottom_left_right_edges(top_depth=top_depth, bottom_depth=bottom_depth, tolerance=tolerance, min_separation=min_separation)
        if top_edges is None or bottom_edges is None or left_edges is None or right_edges is None:
            return None, None, None, None
        
        top_vertices = np.unique(top_edges.flatten())
        bottom_vertices = np.unique(bottom_edges.flatten())
        left_vertices = np.unique(left_edges.flatten())
        right_vertices = np.unique(right_edges.flatten())

        top_tris = np.where(np.isin(self.triangles, top_vertices).any(axis=1))[0]
        bottom_tris = np.where(np.isin(self.triangles, bottom_vertices).any(axis=1))[0]
        left_tris = np.where(np.isin(self.triangles, left_vertices).any(axis=1))[0]
        right_tris = np.where(np.isin(self.triangles, right_vertices).any(axis=1))[0]

        return top_tris, bottom_tris, left_tris, right_tris
    
    def check_intersection(self, other_mesh):
        """
        Checks if this mesh intersects with another mesh.

        Args:
            other_mesh (FaultMesh): Another FaultMesh instance to check for intersection.
        Returns:
            bool: True if the meshes intersect, False otherwise.
        """
        assert isinstance(other_mesh, FaultMesh), "other_mesh must be an instance of FaultMesh."
        assert self.mesh is not None, "This mesh is not defined."
        assert other_mesh.mesh is not None, "Other mesh is not defined."

        pv_self = pv.PolyData(self.vertices, np.hstack([np.full((len(self.triangles), 1), 3), self.triangles]).astype(np.int64))
        pv_other = pv.PolyData(other_mesh.vertices, np.hstack([np.full((len(other_mesh.triangles), 1), 3), other_mesh.triangles]).astype(np.int64))

        intersection = pv_self.boolean_intersection(pv_other, tolerance=1.0e-5)

        return intersection.n_points > 0