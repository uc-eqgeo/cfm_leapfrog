import meshio
import numpy as np
import os
from fault_mesh.utilities.meshing import fit_plane_to_points, get_strike_dip_from_normal
from itertools import combinations
import pyvista as pv
from scipy.spatial import KDTree
from shapely.geometry import box
from scipy.spatial import ConvexHull
from numba import njit

@njit
def cross_3d(a, b):
    return np.array([
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ])

@njit
def triangle_normal(v0, v1, v2):
    """Calculate the normal vector of a triangle defined by vertices v0, v1, v2."""
    edge1 = v1 - v0
    edge2 = v2 - v0
    normal = cross_3d(edge1, edge2)
    norm = np.linalg.norm(normal)
    if norm == 0:
        return normal
    return normal / norm

@njit
def check_if_triangle_crosses_other_triangle(tri1, tri2):
    """Check if two triangles in 3D space intersect."""
    # Using Möller–Trumbore intersection algorithm
    v0, v1, v2 = tri1
    u0, u1, u2 = tri2

    # Compute triangle edges
    e1 = v1 - v0
    e2 = v2 - v0
    f1 = u1 - u0
    f2 = u2 - u0

    # Compute normals
    n1 = triangle_normal(v0, v1, v2)
    n2 = triangle_normal(u0, u1, u2)

    # Check if triangles are coplanar
    if np.abs(np.dot(n1, n2)) == 1.0:
        return False  # Coplanar triangles are not considered intersecting here

    # Compute direction of the line of intersection
    d = cross_3d(n1, n2)
    if np.linalg.norm(d) == 0:
        return False  # Parallel planes

    # Project triangles onto the line of intersection and check for overlap
    def project_onto_line(tri, line_dir):
        projections = [np.dot(vertex, line_dir) for vertex in tri]
        return min(projections), max(projections)

    min1, max1 = project_onto_line(tri1, d)
    min2, max2 = project_onto_line(tri2, d)

    # Check for overlap in projections
    return not (max1 < min2 or max2 < min1)

def find_furthest_points_hull(points):
    """
    Find the two points that are furthest apart using convex hull.
    This is faster for large point clouds as only hull vertices need to be checked.
    
    Parameters:
    -----------
    points : np.ndarray
        Array of shape (n, 3) containing 3D points
        
    Returns:
    --------
    tuple
        (index1, index2, max_distance) - indices of furthest points and their distance
    """
    # Get convex hull vertices
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    
    # Compute distances only between hull vertices
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(hull_points))
    
    # Find the maximum distance
    max_idx = np.unravel_index(np.argmax(distances), distances.shape)
    max_distance = distances[max_idx]
    
    # Map back to original indices
    original_idx1 = hull.vertices[max_idx[0]]
    original_idx2 = hull.vertices[max_idx[1]]
    
    return original_idx1, original_idx2, max_distance


def find_furthest_points_brute(points):
    """
    Find the two points that are furthest apart using brute force.
    
    Parameters:
    -----------
    points : np.ndarray
        Array of shape (n, 3) containing 3D points
        
    Returns:
    --------
    tuple
        (index1, index2, max_distance) - indices of furthest points and their distance
    """
    from scipy.spatial.distance import pdist, squareform
    
    # Compute all pairwise distances
    distances = squareform(pdist(points))
    
    # Find the maximum distance
    max_idx = np.unravel_index(np.argmax(distances), distances.shape)
    max_distance = distances[max_idx]
    
    return max_idx[0], max_idx[1], max_distance


def rectangles_intersect(rect1, rect2):
    """
    Check if two axis-aligned rectangles intersect.
    
    Parameters:
    -----------
    rect1, rect2 : tuple or list
        Rectangle defined as (xmin, ymin, xmax, ymax)
        
    Returns:
    --------
    bool
        True if rectangles intersect, False otherwise
    """
    xmin1, ymin1, xmax1, ymax1 = rect1
    xmin2, ymin2, xmax2, ymax2 = rect2
    
    # Rectangles do NOT intersect if any of these conditions is true:
    # - rect1 is completely to the left of rect2
    # - rect1 is completely to the right of rect2
    # - rect1 is completely above rect2
    # - rect1 is completely below rect2
    
    return not (xmax1 < xmin2 or xmax2 < xmin1 or ymax1 < ymin2 or ymax2 < ymin1)


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
        self._tri_centres = None
        self._tri_centre_dict = None
        self.pv_triangles = None
        self.pv_tri_dict = None
        self.pv_vertices = None
        self.pv_vertex_dict = None
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
        self._tree = None
        self.xmin = None
        self.xmax = None
        self.ymin = None    
        self.ymax = None


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
        pv_mesh = pv.from_meshio(mesh).extract_surface()
        fault_mesh.pv_triangles = pv_mesh.faces.reshape(-1, 4)[:, 1:]
        fault_mesh.vertices = mesh.points
        fault_mesh.pv_vertices = pv_mesh.points
        fault_mesh.triangles = mesh.cells_dict["triangle"]
        fault_mesh.tri_dict = {i: tri for i, tri in enumerate(fault_mesh.triangles)}
        fault_mesh.vertex_dict = {i: vertex for i, vertex in enumerate(fault_mesh.vertices)}
        fault_mesh.pv_tri_dict = {i: tri for i, tri in enumerate(fault_mesh.pv_triangles)}
        fault_mesh.pv_vertex_dict = {i: vertex for i, vertex in enumerate(pv_mesh.points)}
        fault_mesh.name = os.path.splitext(os.path.basename(file_path))[0]
        fault_mesh.xmin = np.min(fault_mesh.vertices[:, 0])
        fault_mesh.xmax = np.max(fault_mesh.vertices[:, 0])
        fault_mesh.ymin = np.min(fault_mesh.vertices[:, 1])
        fault_mesh.ymax = np.max(fault_mesh.vertices[:, 1])
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
            self._resolved_vertex_x = np.dot(self.pv_vertices - self.plane_origin, self.plane_along_strike)
        return self._resolved_vertex_x

    @property
    def resolved_vertex_y(self):
        if self._resolved_vertex_y is None:
            assert self.vertices is not None, "Mesh vertices are not defined."
            self._resolved_vertex_y = np.dot(self.pv_vertices - self.plane_origin, self.plane_down_dip)
        return self._resolved_vertex_y
    
    @property
    def resolved_vertex_z(self):
        if self._resolved_vertex_z is None:
            assert self.vertices is not None, "Mesh vertices are not defined."
            self._resolved_vertex_z = np.dot(self.vertices - self.plane_origin, self.plane_normal)
        return self._resolved_vertex_z

    def build_kdtree(self):
        if self.vertices is not None:
            self._tree = KDTree(self.vertices)
        else:
            self._tree = None

    @property
    def tree(self):
        if self._tree is None:
            self.build_kdtree()
        return self._tree

    def find_all_outside_edges(self):
        triangle_edges = np.array([np.sort(np.array(list(combinations(tri, 2))), axis=1) for tri in self.pv_triangles])
        unique_edges, edge_counts = np.unique(triangle_edges.reshape(-1, 2), axis=0, return_counts=True)
        outside_edges = unique_edges[edge_counts == 1]
        return outside_edges

    def find_all_outside_vertex_indices(self):
        outside_edges = self.find_all_outside_edges()
        outside_vertex_indices = np.unique(outside_edges.reshape(-1, 2))
        return outside_vertex_indices

    def find_all_outside_vertices(self):
        outside_vertex_indices = self.find_all_outside_vertex_indices()
        return self.pv_vertices[outside_vertex_indices]
    
    def find_vertex_indices_constant_depth(self, bottom_value: float, tolerance=10.0):
        outside_vertex_indices = self.find_all_outside_vertex_indices()
        bottom_vertex_indices = np.where(np.abs(self.pv_vertices[:, -1][outside_vertex_indices] - bottom_value) <= tolerance)[0]
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
        for i, tri in self.pv_tri_dict.items():
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
            all_vertex_indices = self.find_all_outside_vertex_indices()
            bottom_vertex_indices = all_vertex_indices[~np.isin(all_vertex_indices, top_vertex_indices)]
            return top_vertex_indices, bottom_vertex_indices, None, None
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
        # if top_vertex_indices is None or bottom_vertex_indices is None or left_vertex_indices is None or right_vertex_indices is None:
        #     return None, None, None, None
        all_edges = self.find_all_outside_edges()
        top_edges = []
        bottom_edges = []
        if left_vertex_indices is None or right_vertex_indices is None:
            left_edges = None
            right_edges = None 
        else:
            left_edges = []
            right_edges = []

        for edge in all_edges:
            if edge[0] in top_vertex_indices and edge[1] in top_vertex_indices:
                top_edges.append(edge)
            elif edge[0] in bottom_vertex_indices and edge[1] in bottom_vertex_indices:
                bottom_edges.append(edge)
            if left_edges is not None and right_edges is not None:
                if any([vert_i in left_vertex_indices for vert_i in edge]):
                        left_edges.append(edge)
                elif any([vert_i in right_vertex_indices for vert_i in edge]):
                        right_edges.append(edge)

        if left_edges is not None:
            left_edges = np.array(left_edges)
        if right_edges is not None:
            right_edges = np.array(right_edges)

        return np.array(top_edges), np.array(bottom_edges), left_edges, right_edges
    
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
        if top_edges is None or bottom_edges is None:
            return None, None, None, None
        
        top_vertices = np.unique(top_edges.flatten())
        bottom_vertices = np.unique(bottom_edges.flatten())

        top_tris = np.where(np.isin(self.pv_triangles, top_vertices).any(axis=1))[0]
        bottom_tris = np.where(np.isin(self.pv_triangles, bottom_vertices).any(axis=1))[0]

        if left_edges is None or right_edges is None:
            left_tris = None
            right_tris = None
        else:
            left_vertices = np.unique(left_edges.flatten())
            right_vertices = np.unique(right_edges.flatten())
            left_tris = np.where(np.isin(self.pv_triangles, left_vertices).any(axis=1))[0]
            right_tris = np.where(np.isin(self.pv_triangles, right_vertices).any(axis=1))[0]

        return top_tris, bottom_tris, left_tris, right_tris
    
    def check_intersection2d(self, other_mesh):
        """
        Checks if this mesh intersects with another mesh in 2D (ignoring depth).

        Args:
            other_mesh (FaultMesh): Another FaultMesh instance to check for intersection.
        Returns:
            bool: True if the meshes intersect in 2D, False otherwise.
        """
        assert isinstance(other_mesh, FaultMesh), "other_mesh must be an instance of FaultMesh."
        assert self.mesh is not None, "This mesh is not defined."
        assert other_mesh.mesh is not None, "Other mesh is not defined."

        # Quick bounding box check
        return rectangles_intersect(
            (self.xmin, self.ymin, self.xmax, self.ymax),
            (other_mesh.xmin, other_mesh.ymin, other_mesh.xmax, other_mesh.ymax)
        )
    
    def check_intersection(self, other_mesh, min_kd_distance=3.e3):
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

        if not self.check_intersection2d(other_mesh):
            return False

        # Quick check using KDTree to see if any vertices are within min_kd_distance
        distances, _ = self.tree.query(other_mesh.vertices, k=1)
        if np.all(distances > min_kd_distance):
            return False
        
        else:
            return True
        
        # # Use PyVista for a more precise intersection check

        # pv_self = pv.from_meshio(self.mesh).extract_surface()
        # pv_other = pv.from_meshio(other_mesh.mesh).extract_surface()

        # clipped = pv_self.clip_surface(pv_other, invert=False)
        # if clipped.n_cells == pv_self.n_cells:
        #     return False
        # else:
        #     return True
        
    def get_intersecting_triangles(self, other_mesh, max_distance=3.e3, max_implicit_distance: float = None):
        """
        Gets the triangles from this mesh that intersect with another mesh.

        Args:
            other_mesh (FaultMesh): Another FaultMesh instance to check for intersection.

        Returns:
            np.ndarray: Array of triangle indices from this mesh that intersect with the other mesh.
        """
        assert isinstance(other_mesh, FaultMesh), "other_mesh must be an instance of FaultMesh."
        assert self.mesh is not None, "This mesh is not defined."
        assert other_mesh.mesh is not None, "Other mesh is not defined."

        if not self.check_intersection(other_mesh):
            return np.array([], dtype=int)

        # Use pyvista for intersection check
        pv_self = pv.from_meshio(self.mesh).extract_surface()
        pv_other = pv.from_meshio(other_mesh.mesh).extract_surface()

        tri_centres_self = np.mean(pv_self.points[pv_self.faces.reshape(-1, 4)[:, 1:]], axis=1)
        tri_centres_other = np.mean(pv_other.points[pv_other.faces.reshape(-1, 4)[:, 1:]], axis=1)

        kdtree = KDTree(tri_centres_other)
        distances, _ = kdtree.query(tri_centres_self, k=1)
        close_tri_indices = np.where(distances < max_distance)[0]


        implicit_distances = np.array(pv_self.compute_implicit_distance(pv_other)["implicit_distance"])
        triangle_distances = implicit_distances[self.pv_triangles]
        if max_implicit_distance is not None:
            far_tri_indices = np.where(np.abs(triangle_distances) < max_implicit_distance)[0]
            close_tri_indices = np.intersect1d(close_tri_indices, far_tri_indices)
        else:
            all_positive = np.all(triangle_distances >= 0, axis=1)
            all_negative = np.all(triangle_distances < 0, axis=1)
            to_keep = ~(all_positive | all_negative)
            where_true = np.where(to_keep)[0]
            close_tri_indices = np.intersect1d(where_true, close_tri_indices)

        return close_tri_indices
    
    def generate_cutting_mesh(self, other_mesh, max_distance=3.e3, max_implicit_distance: float = None):
        intersecting_triangles = self.get_intersecting_triangles(other_mesh, max_distance=max_distance, max_implicit_distance=max_implicit_distance)
        if len(intersecting_triangles) == 0:
            return None
        
        cutting_vertices = []
        cutting_triangles = []
        vertex_map = {}
        for tri_index in intersecting_triangles:
            tri_vertices = self.pv_triangles[tri_index]
            new_tri_vertices = []
            for v in tri_vertices:
                if v not in vertex_map:
                    vertex_map[v] = len(cutting_vertices)
                    cutting_vertices.append(self.pv_vertex_dict[v])
                new_tri_vertices.append(vertex_map[v])
            cutting_triangles.append(new_tri_vertices)
        
        cutting_mesh = meshio.Mesh(
            points=np.array(cutting_vertices),
            cells=[("triangle", np.array(cutting_triangles))]
        )
        return cutting_mesh
    

    def decide_whether_to_cut(self, other_mesh, threshold=0.9, min_distance: float = 10.e3, higher_meshes: list = None, higher_mesh_tolerance: float = 1.e3):
        intersecting_triangles = self.get_intersecting_triangles(other_mesh)
        if len(intersecting_triangles) == 0:
            return False
        
        # Check whether faults are already cut by higher priority faults
        if higher_meshes is not None:
            for higher_mesh in higher_meshes:
                self_check_intersection = self.check_intersection(higher_mesh)
                other_check_intersection = other_mesh.check_intersection(higher_mesh)
                if self_check_intersection and other_check_intersection:
                    self_implicit_distances = np.array(pv.from_meshio(self.mesh).extract_surface().compute_implicit_distance(pv.from_meshio(higher_mesh.mesh).extract_surface())["implicit_distance"])
                    other_implicit_distances = np.array(pv.from_meshio(other_mesh.mesh).extract_surface().compute_implicit_distance(pv.from_meshio(higher_mesh.mesh).extract_surface())["implicit_distance"])

                    self_mainly_positive = np.sum(self_implicit_distances > 0) > 0.7 * len(self_implicit_distances)
                    other_mainly_positive = np.sum(other_implicit_distances > 0) > 0.7 * len(other_implicit_distances)

                    if self_mainly_positive != other_mainly_positive:
                        if self_mainly_positive:
                            self_gt_threshold = np.sum(self_implicit_distances < -higher_mesh_tolerance)
                            other_gt_threshold = np.sum(other_implicit_distances >= higher_mesh_tolerance)
                        else:
                            self_gt_threshold = np.sum(self_implicit_distances >= higher_mesh_tolerance)
                            other_gt_threshold = np.sum(other_implicit_distances < -higher_mesh_tolerance)
                        if not any([self_gt_threshold, other_gt_threshold]):
                            print(f"Deciding NOT to cut {self.name} against {other_mesh.name} due to higher priority fault {higher_mesh.name}.")
                            return False
        

        
        top_tris, bottom_tris, left_tris, right_tris = self.find_top_bottom_left_right_triangles()
        edge_conditions = {
            'top': len(np.intersect1d(intersecting_triangles, top_tris)) > 0,
            'bottom': len(np.intersect1d(intersecting_triangles, bottom_tris)) > 0,
        }
        if left_tris is not None and right_tris is not None:
            edge_conditions['left'] = len(np.intersect1d(intersecting_triangles, left_tris)) > 0
            edge_conditions['right'] = len(np.intersect1d(intersecting_triangles, right_tris)) > 0
        if sum(edge_conditions.values()) >= 2:
            return True
        
        
        top_vertex_indices, bottom_vertex_indices, left_vertex_indices, right_vertex_indices = self.find_top_bottom_left_right_vertices()
        top_vertex_points, bottom_vertex_points = self.pv_vertices[top_vertex_indices], self.pv_vertices[bottom_vertex_indices]
        if left_vertex_indices is None or right_vertex_indices is None:
            left_vertex_points, right_vertex_points = None, None
        else:
            left_vertex_points, right_vertex_points = self.pv_vertices[left_vertex_indices], self.pv_vertices[right_vertex_indices]
        vertex_point_dict = {
            'top': top_vertex_points, 
            'bottom': bottom_vertex_points,
            'left': left_vertex_points,
            'right': right_vertex_points
        }
        # Find centres of intersecting triangles
        centres = np.mean(self.pv_vertices[self.pv_triangles[intersecting_triangles]], axis=1)
        
        # Find centre points that are furthest apart
        if len(centres) < 2:
            return False
        
        # Use convex hull method for efficiency (or brute force for small arrays)
        if len(centres) > 100:
            idx1, idx2, max_distance = find_furthest_points_hull(centres)
        else:
            idx1, idx2, max_distance = find_furthest_points_brute(centres)
        
        centre1, centre2 = centres[idx1], centres[idx2]

        centre1_nearest_edge = {
            'top': np.min(np.linalg.norm(top_vertex_points - centre1, axis=1)),
            'bottom': np.min(np.linalg.norm(bottom_vertex_points - centre1, axis=1)),
        }
        if left_vertex_points is not None and right_vertex_points is not None:
            centre1_nearest_edge['left'] = np.min(np.linalg.norm(left_vertex_points - centre1, axis=1))
            centre1_nearest_edge['right'] = np.min(np.linalg.norm(right_vertex_points - centre1, axis=1))
        centre2_nearest_edge = {
            'top': np.min(np.linalg.norm(top_vertex_points - centre2, axis=1)),
            'bottom': np.min(np.linalg.norm(bottom_vertex_points - centre2, axis=1)),
        }
        if left_vertex_points is not None and right_vertex_points is not None:
            centre2_nearest_edge['left'] = np.min(np.linalg.norm(left_vertex_points - centre2, axis=1))
            centre2_nearest_edge['right'] = np.min(np.linalg.norm(right_vertex_points - centre2, axis=1))
        centre1_closest_edge = min(centre1_nearest_edge, key=centre1_nearest_edge.get)
        centre2_closest_edge = min(centre2_nearest_edge, key=centre2_nearest_edge.get)

        centre1_closest_edge_point = vertex_point_dict[centre1_closest_edge][np.argmin(np.linalg.norm(vertex_point_dict[centre1_closest_edge][:, :2] - centre1[:2], axis=1))]
        centre2_closest_edge_point = vertex_point_dict[centre2_closest_edge][np.argmin(np.linalg.norm(vertex_point_dict[centre2_closest_edge][:, :2] - centre2[:2], axis=1))]
        edge_dist = np.linalg.norm(centre1_closest_edge_point - centre2_closest_edge_point)

        print(f"Furthest centres are {max_distance/1.e3:.1f} km apart, closest to edges {centre1_closest_edge} and {centre2_closest_edge}.")
        if centre1_closest_edge != centre2_closest_edge:
            if (max_distance / edge_dist) > threshold:
                print(f"Deciding to cut {self.name} against {other_mesh.name}: max_distance / edge_dist = {max_distance/edge_dist:.2f} > {threshold}")
                return True
            else:
                print(f"Deciding NOT to cut {self.name} against {other_mesh.name}: max_distance / edge_dist = {max_distance/edge_dist:.2f} <= {threshold}")
                return False
        else:
            _, _, edge_length = find_furthest_points_brute(vertex_point_dict[centre1_closest_edge])
            if edge_length < min_distance:
                if (max_distance / edge_length) > threshold:
                    print(f"Deciding to cut {self.name} against {other_mesh.name}: max_distance / edge_length = {max_distance/edge_length:.2f} > {threshold}")
                    return True
                else:
                    print(f"Deciding NOT to cut {self.name} against {other_mesh.name}: max_distance / edge_length = {max_distance/edge_length:.2f} <= {threshold}")
                    return False
            else:
                if all([centre1_nearest_edge[centre1_closest_edge] < 3.e3, centre2_nearest_edge[centre2_closest_edge] < 3.e3]):
                    print(f"Both centres are very close to the same edge: {centre1_closest_edge}")
                    if edge_dist > min_distance:
                        print(f"Deciding to cut {self.name} against {other_mesh.name}: edge_length = {edge_length:.2f} > {min_distance}")
                        return True
                    
                    elif (edge_dist/edge_length) > 0.5:
                        print(f"Deciding to cut {self.name} against {other_mesh.name}: edge_dist / edge_length = {edge_dist/edge_length:.2f} > 0.5")
                        return True

                    else:
                        print(f"Deciding NOT to cut {self.name} against {other_mesh.name}: edge_length = {edge_length:.2f} <= {min_distance}")
                        return False
                else:
                    if any([centre1_nearest_edge[centre1_closest_edge] < 3.e3, centre2_nearest_edge[centre2_closest_edge] < 3.e3]):
                        print(f"At least one centre is not very close to the shared edge: {centre1_closest_edge}")
                        if centre1_nearest_edge[centre1_closest_edge] < 3.e3:
                            centre2_second_closest_edge = sorted(centre2_nearest_edge.items(), key=lambda item: item[1])[1][0]
                            centre2_second_closest_edge_point = vertex_point_dict[centre2_second_closest_edge][np.argmin(np.linalg.norm(vertex_point_dict[centre2_second_closest_edge][:, :2] - centre2[:2], axis=1))]
                            edge_dist_second = np.linalg.norm(centre1_closest_edge_point - centre2_second_closest_edge_point)
                            print(f"Centre 2 second closest edge: {centre2_second_closest_edge} at distance {centre2_nearest_edge[centre2_second_closest_edge]}")
                            if (max_distance / edge_dist_second) > threshold:
                                print(f"Deciding to cut {self.name} against {other_mesh.name}: max_distance / edge_dist_second = {max_distance/edge_dist_second:.2f} > {threshold}")
                                return True
                            else:
                                print(f"Deciding NOT to cut {self.name} against {other_mesh.name}: max_distance / edge_dist_second = {max_distance/edge_dist_second:.2f} <= {threshold}")
                                return False
                        else:
                            centre1_second_closest_edge = sorted(centre1_nearest_edge.items(), key=lambda item: item[1])[1][0]
                            centre1_second_closest_edge_point = vertex_point_dict[centre1_second_closest_edge][np.argmin(np.linalg.norm(vertex_point_dict[centre1_second_closest_edge][:, :2] - centre1[:2], axis=1))]
                            edge_dist_second = np.linalg.norm(centre1_second_closest_edge_point - centre2_closest_edge_point)
                            print(f"Centre 1 second closest edge: {centre1_second_closest_edge} at distance {edge_dist_second}")
                            if (max_distance / edge_dist_second) > threshold:
                                print(f"Deciding to cut {self.name} against {other_mesh.name}: max_distance / edge_dist_second = {max_distance/edge_dist_second:.2f} > {threshold}")
                                return True
                            else:
                                print(f"Deciding NOT to cut {self.name} against {other_mesh.name}: max_distance / edge_dist_second = {max_distance/edge_dist_second:.2f} <= {threshold}")
                                return False
                    else:
                        print(f"Neither centre is very close to any edge: NOT cutting.")
                        return False
                
    def cut_mesh(self, other_mesh, fault_trace: np.ndarray, surface_tolerance = 1.0, cutting_fragment: meshio.Mesh = None):
        """
        Cuts this mesh using another mesh.

        Args:
            other_mesh (FaultMesh): Another FaultMesh instance to cut this mesh with.

        Returns:
            FaultMesh: A new FaultMesh instance representing the cut mesh.
        """
        assert isinstance(other_mesh, FaultMesh), "other_mesh must be an instance of FaultMesh or meshio.Mesh."
        assert self.mesh is not None, "This mesh is not defined."
        assert other_mesh.mesh is not None, "Other mesh is not defined."
        assert fault_trace.shape[1] == 3, "fault_trace must be an array of shape (n, 3)."

        if not self.decide_whether_to_cut(other_mesh):
            print(f"Decided not to cut {self.name} with {other_mesh.name}.")
            return self

        pv_self = pv.from_meshio(self.mesh).extract_surface()
        if cutting_fragment is not None:
            pv_other = pv.from_meshio(cutting_fragment).extract_surface()
        else:
            pv_other = pv.from_meshio(other_mesh.mesh).extract_surface()

        clipped1 = pv_self.clip_surface(pv_other, invert=False)
        clipped1_surface_points = clipped1.points[clipped1.points[:, 2] >= -surface_tolerance]
        clipped2 = pv_self.clip_surface(pv_other, invert=True)
        clipped2_surface_points = clipped2.points[clipped2.points[:, 2] >= -surface_tolerance]

        if clipped1_surface_points.shape[0] == 0:
            print(f"Clipped1 resulted in no surface points for {self.name} after cutting with {other_mesh.name}. Returning clipped2.")
            clipped = clipped2
        if clipped2_surface_points.shape[0] == 0:
            print(f"Clipped2 resulted in no surface points for {self.name} after cutting with {other_mesh.name}. Returning clipped1.")
            clipped = clipped1

        else:
            trace_kdtree = KDTree(fault_trace)
            clipped_1_dist, _ = trace_kdtree.query(clipped1_surface_points)
            clipped_2_dist, _ = trace_kdtree.query(clipped2_surface_points)

            dist1 = np.mean(clipped_1_dist)
            dist2 = np.mean(clipped_2_dist)

            if dist1 < dist2:
                clipped = clipped1
            else:
                clipped = clipped2

        if clipped.n_points < 3 or clipped.n_cells < 1:
            print(f"Cutting resulted in too few points or cells for {self.name} after cutting with {other_mesh.name}. Returning original mesh.")
            return self

        new_fault_mesh = FaultMesh()
        new_fault_mesh.mesh = meshio.Mesh(points=clipped.points, cells={"triangle": clipped.faces.reshape(-1, 4)[:, 1:]})
        new_fault_mesh.name = f"{self.name}_cut_by_{other_mesh.name}"
        new_fault_mesh.vertices = new_fault_mesh.mesh.points
        new_fault_mesh.triangles = new_fault_mesh.mesh.cells_dict["triangle"]
        new_fault_mesh.tri_dict = {i: tri for i, tri in enumerate(new_fault_mesh.triangles)}
        new_fault_mesh.pv_triangles = clipped.faces.reshape(-1, 4)[:, 1:]
        new_fault_mesh.pv_tri_dict = {i: tri for i, tri in enumerate(new_fault_mesh.pv_triangles)}
        new_fault_mesh.vertex_dict = {i: vertex for i, vertex in enumerate(new_fault_mesh.vertices)}
        new_fault_mesh.pv_vertices = clipped.points
        new_fault_mesh.pv_vertex_dict = {i: vertex for i, vertex in enumerate(new_fault_mesh.pv_vertices)}
        new_fault_mesh.name = f"{self.name}"
        new_fault_mesh.xmin = np.min(new_fault_mesh.vertices[:, 0])
        new_fault_mesh.xmax = np.max(new_fault_mesh.vertices[:, 0])
        new_fault_mesh.ymin = np.min(new_fault_mesh.vertices[:, 1])
        new_fault_mesh.ymax = np.max(new_fault_mesh.vertices[:, 1])

        return new_fault_mesh
                


