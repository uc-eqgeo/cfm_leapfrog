# Mesh cutting logic

When two fault meshes intersect or nearly intersect, one mesh is cut by the other so that their geometries are consistent — the lower-priority fault terminates at the surface of the higher-priority fault. This page describes how the code decides whether a cut is needed and, once it is, which half of the cut mesh to keep.

## Deciding whether to cut

The entry point is `FaultMesh.decide_whether_to_cut(self, other_mesh, ...)`.

### Step 1: check for geometric intersection

The method first calls `get_intersecting_triangles(other_mesh)` to find triangles in `self` that intersect or come within a threshold distance of `other_mesh`. If no such triangles exist the meshes do not interact and the method returns `False` immediately — no cut is performed.

### Step 2: override by higher-priority faults

If a list of `higher_meshes` is supplied, the code checks whether a higher-priority fault has already separated `self` and `other_mesh` onto opposite sides of itself. For each higher-priority mesh:

1. Both `self` and `other_mesh` must intersect the higher-priority mesh (otherwise the check is skipped for that pair).
2. The signed implicit distances from every vertex of `self` and `other_mesh` to the higher-priority mesh surface are computed.
3. If the two faults lie predominantly on *opposite* sides of the higher-priority mesh (one is mainly positive, the other mainly negative) **and** neither fault has significant penetration beyond the higher-priority surface (i.e. neither has many vertices more than `higher_mesh_tolerance` on the wrong side), then the higher-priority fault is already acting as a separator. In this situation, cutting `self` against `other_mesh` is unnecessary and the method returns `False`.

### Step 3: edge-intersection test (always active)

The intersecting triangles are compared against the four edge groups of `self` (top, bottom, left, right — identified by depth and by the largest spatial gap in the side vertices). If the intersection zone touches **two or more distinct edges** the intersection genuinely crosses the mesh, and the method returns `True`.

This test catches the common case where one fault clearly cuts through another from one edge to another (e.g. top to bottom, or left to right).

### Step 4: fancy cutting heuristics (optional)

When `fancy_cutting=True` and the simple edge test did not fire, additional heuristics are applied.

1. **Find the two most widely separated intersection-triangle centres.** These are the centres of the intersecting triangles that are furthest apart from each other (`max_distance`). For large intersection zones the convex-hull method is used for efficiency; for small zones brute-force pairwise distances are used.

2. **Identify the closest mesh edge for each centre.** Each of the two extreme centres is matched to its nearest mesh edge (top, bottom, left, or right) by minimum Euclidean distance to the edge vertices.

3. **Decision tree:**

   - **Centres are closest to different edges** — the intersection zone spans two different edges, which strongly suggests a genuine through-cutting geometry. The ratio `max_distance / edge_dist` is compared to `threshold` (default 0.9), where `edge_dist` is the distance between the closest edge points to each centre. If the ratio exceeds the threshold, cut; otherwise do not.

   - **Centres are closest to the same edge, and that edge is short** (`edge_length < min_distance`) — the intersection is near a short edge. The ratio `max_distance / edge_length` is compared to `threshold`. If it exceeds the threshold, cut.

   - **Centres are closest to the same edge, and that edge is long:**
     - If both centres are very close to that edge (within 3 km), the intersection zone hugs a single long edge. Cut if `edge_dist > min_distance` or `edge_dist / edge_length > 0.5`; otherwise do not cut.
     - If only one centre is very close to the edge, use the *second*-closest edge for that centre instead and re-compute `edge_dist`. Apply the threshold test on `max_distance / edge_dist_second`.
     - If neither centre is close to any edge, do not cut.

   If `fancy_cutting=False` (the default), the method returns `False` after the edge test.

---

## Choosing which side to keep after cutting

Once a cut is decided, `FaultMesh.cut_mesh` (or `cut_mesh_pv`) clips `self` against `other_mesh` in both directions using PyVista's `clip_surface`:

- `clipped1` — the portion of `self` that lies on the *outside* of `other_mesh` (`invert=False`)
- `clipped2` — the portion of `self` that lies on the *inside* of `other_mesh` (`invert=True`)

### Degenerate cases

If either clipped half has no surface points (points with z ≥ −`surface_tolerance`, i.e. near the Earth's surface), it is discarded immediately and the other half is kept. This handles cases where the clip removes essentially all of one side.

### General case: proximity to the fault trace

In the normal case both halves have surface points. A KD-tree is built from the supplied `fault_trace` — the 3-D polyline of the fault's surface expression.

**`cut_mesh`** counts, for each half, how many surface points lie within `surface_max_dist_tolerance` of the trace. The half with *more* nearby points is kept. This reflects the expectation that the correct side of the cut mesh should have its surface expression close to the mapped fault trace.

**`cut_mesh_pv`** uses a slightly different criterion: it keeps the half whose surface points have the *lower mean distance* to the fault trace, which is more robust when the two halves have different numbers of surface points.

### Safety check

If the retained half has fewer than 3 points or 1 cell (i.e. the clip produced an essentially empty mesh), the original uncut mesh is returned unchanged with a warning.
