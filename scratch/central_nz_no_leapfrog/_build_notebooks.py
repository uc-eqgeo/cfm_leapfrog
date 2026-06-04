"""Generate the combined central-NZ fault-mesh pipeline notebook in English
and Japanese. Run once; produces two .ipynb files, then can be deleted."""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

# Each entry: ("md"|"code", english_text, japanese_text).
# Code cells are identical apart from inline comments (translated).
CELLS = []
def md(en, ja):  CELLS.append(("md", en, ja))
def code(en, ja): CELLS.append(("code", en, ja))

# ---------------------------------------------------------------- title
md(
"""# Central New Zealand fault-mesh pipeline

This notebook builds triangulated fault-surface meshes for the central New
Zealand Community Fault Model, from raw fault traces all the way to clean,
uniformly-sized triangle meshes. It combines what used to be three separate
scripts into a single, ordered workflow:

| Stage | What it does | Outputs |
|-------|--------------|---------|
| **1. Build surfaces** | Turn each (multi-segment) fault into a 3-D triangulated surface from depth contours | `test_objs/`, `test_vtks/`, `test_contours/` |
| **2. Cut surfaces** | Trim each surface where higher-priority faults cross it, and against the base depth surface | `test_final_meshes/*_cut.obj` |
| **3. Remesh** | Optimise every cut surface into uniform, near-equilateral triangles with **MMG** | `test_remeshed_vtks/`, `merged_mesh_remeshed.vtk` |

**Run the cells in order, top to bottom.** The stages are sequential: each one
consumes the files written by the previous one, so a single linear run produces
everything. You can also re-run a later stage on its own (e.g. re-cut without
re-building) because each stage reads its inputs back from disk.

> **Prerequisites**: run this notebook from its own folder
> (`scratch/central_nz_no_leapfrog/`) inside the `leapfrog-fault-models` conda
> environment. All paths below are relative to that folder.
""",
"""# 中央ニュージーランド断層メッシュ・パイプライン

このノートブックは、ニュージーランド中央部のコミュニティ断層モデル（Community
Fault Model）について、断層トレースから出発して、最終的に均一なサイズの三角形
メッシュに整えるまでの、断層面の三角形メッシュ生成を一貫して行います。これまで
3つに分かれていたスクリプトを、順序立てた1つのワークフローにまとめたものです：

| ステージ | 内容 | 出力 |
|-------|--------------|---------|
| **1. 面の生成** | 各（複数セグメントの）断層を、深度コンターから3次元の三角形面に変換 | `test_objs/`, `test_vtks/`, `test_contours/` |
| **2. 面のカット** | 上位の断層が横切る箇所、および基底の深度面に沿って各面をトリミング | `test_final_meshes/*_cut.obj` |
| **3. リメッシュ** | カット済みの各面を **MMG** で均一・ほぼ正三角形のメッシュに最適化 | `test_remeshed_vtks/`, `merged_mesh_remeshed.vtk` |

**セルは上から順番に実行してください。** 各ステージは前のステージが書き出した
ファイルを入力として使う逐次処理です。したがって上から一通り実行すれば、すべて
の成果物が生成されます。各ステージは入力をディスクから読み直すため、後段だけを
個別に再実行することもできます（例：再生成せずに再カットのみ行う）。

> **前提条件**：このノートブックは、`leapfrog-fault-models` conda 環境で、
> ノートブック自身のフォルダ（`scratch/central_nz_no_leapfrog/`）から実行して
> ください。以下のパスはすべてそのフォルダからの相対パスです。
""")

# ---------------------------------------------------------------- imports
md(
"""## Setup

### Imports""",
"""## セットアップ

### インポート""")

code(
"""import os
from pathlib import Path

import numpy as np
import meshio
import pyvista as pv

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.mesh import FaultMesh
from fault_mesh.io.array_operations import read_raster""",
"""import os
from pathlib import Path

import numpy as np
import meshio
import pyvista as pv

from fault_mesh.faults.leapfrog import LeapfrogMultiFault
from fault_mesh.faults.mesh import FaultMesh
from fault_mesh.io.array_operations import read_raster""")

# ---------------------------------------------------------------- config
md(
"""### Configuration

Every input path and tunable parameter lives in this one cell, so a collaborator
can adapt the pipeline to a new dataset by editing here and nowhere else.

The `*_edited.csv` files are **human-curated** inputs (fault groupings, cutting
order, trace extensions) — they are reviewed and edited by hand, not generated
automatically. See the comments by each path.""",
"""### 設定

すべての入力パスと調整可能なパラメータは、この1つのセルにまとめてあります。
共同研究者は、新しいデータセットに合わせる際、ここだけを編集すれば対応できます。

`*_edited.csv` ファイルは**人手で精査した**入力です（断層のグループ分け、カット
順序、トレース延長）。自動生成ではなく、手作業でレビュー・編集します。各パスの
コメントを参照してください。""")

code(
"""# --- Input files (relative to this notebook's folder) ---
DOCS = Path("../../docs/tutorials")
FAULT_SHP         = DOCS / "tutorial_gis/central_nz_minimal_data.shp"           # raw fault traces
FAULT_SYSTEMS     = DOCS / "define_connections_data/central_gt1_5_connected_edited.csv"  # which segments form one fault
CUTTING_HIERARCHY = DOCS / "define_connections_data/central_gt1_5_hierarchy_edited.csv"  # which fault cuts which
DEPTH_RASTER      = Path("with_hannu_mods.tif")          # base/Moho-style depth surface (GeoTIFF, z in metres)
ADDITIONAL_CUTS   = Path("additional_cuts.csv")          # extra cut relationships not implied by the hierarchy
# EXCLUDED_CUTS   = Path("excluded_cuts.csv")            # OPTIONAL: pairs that must NEVER be cut (not used here)
TRACE_EXT_EDITED  = Path("trace_extensions_edited.csv")  # curated along-strike trace extensions

# --- Coordinate system & network-building parameters ---
EPSG             = 2193      # NZTM2000; set to None to skip reprojection
DIST_TOLERANCE   = 1000.0    # max gap (m) for two trace segments to be treated as connected

# --- Stage 1: surface meshing ---
DEPTH_CONTOUR_LEVELS = np.arange(0., 32000., 500.)  # depths (m) at which to draw contours
MESH_RESOLUTION      = 500.0                         # target triangle size (m) for the raw surface

# --- Stage 2: cutting ---
MIN_CUT_DISTANCE = 2.0e3     # don't cut faults whose meshes are closer than this (m)
BOTTOM_DEPTH     = -31500.0  # depth (m) used when reasoning about cut extents

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

for d in (OBJ_DIR, VTK_DIR, CONTOUR_DIR, FINAL_DIR, REMESH_DIR):
    d.mkdir(exist_ok=True)""",
"""# --- 入力ファイル（このノートブックのフォルダからの相対パス）---
DOCS = Path("../../docs/tutorials")
FAULT_SHP         = DOCS / "tutorial_gis/central_nz_minimal_data.shp"           # 生の断層トレース
FAULT_SYSTEMS     = DOCS / "define_connections_data/central_gt1_5_connected_edited.csv"  # どのセグメントが1つの断層を成すか
CUTTING_HIERARCHY = DOCS / "define_connections_data/central_gt1_5_hierarchy_edited.csv"  # どの断層がどの断層をカットするか
DEPTH_RASTER      = Path("with_hannu_mods.tif")          # 基底（モホ的）深度面（GeoTIFF、z は m）
ADDITIONAL_CUTS   = Path("additional_cuts.csv")          # 階層からは導かれない追加のカット関係
# EXCLUDED_CUTS   = Path("excluded_cuts.csv")            # 任意：絶対にカットしないペア（ここでは未使用）
TRACE_EXT_EDITED  = Path("trace_extensions_edited.csv")  # 精査済みの走向方向トレース延長

# --- 座標系・ネットワーク構築パラメータ ---
EPSG             = 2193      # NZTM2000。再投影しない場合は None
DIST_TOLERANCE   = 1000.0    # 2つのトレースセグメントを「接続」とみなす最大ギャップ（m）

# --- ステージ1：面メッシュ生成 ---
DEPTH_CONTOUR_LEVELS = np.arange(0., 32000., 500.)  # コンターを描く深度（m）
MESH_RESOLUTION      = 500.0                         # 生の面の目標三角形サイズ（m）

# --- ステージ2：カット ---
MIN_CUT_DISTANCE = 2.0e3     # メッシュ同士がこれより近い断層はカットしない（m）
BOTTOM_DEPTH     = -31500.0  # カット範囲の判定に用いる深度（m）

# --- ステージ3：MMG リメッシュ ---
TARGET_SIZE     = 1000.0     # 均一な目標辺長（m）
HAUSD           = 500.0      # リメッシュ面の元面からの最大許容偏差（m）
RIDGE_ANGLE_DEG = 45.0       # これを超える二面角の折れ目（キンク）を保持

# --- 出力ディレクトリ ---
OBJ_DIR     = Path("test_objs")           # 生の面（OBJ）-- ステージ2の入力
VTK_DIR     = Path("test_vtks")           # 生の面（VTK）-- 確認用
CONTOUR_DIR = Path("test_contours")       # 深度コンター（GeoJSON）
FINAL_DIR   = Path("test_final_meshes")   # カット済みの面 -- ステージ3の入力
REMESH_DIR  = Path("test_remeshed_vtks")  # リメッシュ済みの面
MERGED_VTK  = Path("merged_mesh_remeshed.vtk")

for d in (OBJ_DIR, VTK_DIR, CONTOUR_DIR, FINAL_DIR, REMESH_DIR):
    d.mkdir(exist_ok=True)""")

# ---------------------------------------------------------------- network
md(
"""### Build the fault network

This is the shared foundation for every stage. We:

1. read the raw fault traces from the shapefile and reproject to NZTM;
2. **find connections** — segments whose ends are within `DIST_TOLERANCE` are
   candidates to belong to the same fault;
3. **read fault systems** — the curated CSV says which segments actually form
   one continuous (multi-segment) fault, and we generate those *curated faults*;
4. **read the cutting hierarchy** — the curated CSV that ranks faults so we know,
   later, which fault truncates which.

`fault_data` is built **once** here and reused by all three stages.""",
"""### 断層ネットワークの構築

これは全ステージ共通の土台です。ここでは：

1. シェープファイルから生の断層トレースを読み込み、NZTM に再投影します；
2. **接続を検出** — 端点が `DIST_TOLERANCE` 以内にあるセグメントは、同一断層に
   属する候補となります；
3. **断層システムを読み込み** — 精査済み CSV が、どのセグメントが実際に1本の
   連続した（複数セグメントの）断層を成すかを指定し、それらの*精査済み断層*を
   生成します；
4. **カット階層を読み込み** — 断層に順位を付けた精査済み CSV により、後段でどの
   断層がどの断層を切断するかが分かります。

`fault_data` はここで**一度だけ**構築し、3つのステージすべてで再利用します。""")

code(
"""fault_data = LeapfrogMultiFault.from_nz_cfm_shp(
    str(FAULT_SHP), remove_colons=True, epsg=EPSG,
    smoothing_n=None, dip_multiplier=1.0, exclude_zero=False)
fault_data.segment_distance_tolerance = DIST_TOLERANCE
fault_data.find_connections(verbose=False)

# Curated multi-segment faults, and the order in which faults cut each other.
fault_data.read_fault_systems(str(FAULT_SYSTEMS))
fault_data.generate_curated_faults()
fault_data.read_cutting_hierarchy(str(CUTTING_HIERARCHY))

print(f"{len(fault_data.curated_faults)} curated faults; "
      f"cutting hierarchy has {len(fault_data.cutting_hierarchy)} entries")""",
"""fault_data = LeapfrogMultiFault.from_nz_cfm_shp(
    str(FAULT_SHP), remove_colons=True, epsg=EPSG,
    smoothing_n=None, dip_multiplier=1.0, exclude_zero=False)
fault_data.segment_distance_tolerance = DIST_TOLERANCE
fault_data.find_connections(verbose=False)

# 精査済みの複数セグメント断層と、断層が互いをカットする順序。
fault_data.read_fault_systems(str(FAULT_SYSTEMS))
fault_data.generate_curated_faults()
fault_data.read_cutting_hierarchy(str(CUTTING_HIERARCHY))

print(f"{len(fault_data.curated_faults)} curated faults; "
      f"cutting hierarchy has {len(fault_data.cutting_hierarchy)} entries")""")

# ---------------------------------------------------------------- trace ext
md(
"""### Trace extensions (human-in-the-loop)

Faults often need to be extended a little along strike so that neighbouring
faults actually meet and cut cleanly. This is a **two-step, human-reviewed**
process:

1. `suggest_trace_extensions(...)` writes a CSV (and GeoJSON) of *proposed*
   extensions for you to inspect in GIS.
2. You edit/approve them, save as `trace_extensions_edited.csv`, and
   `read_trace_extensions(...)` + `apply_trace_extensions()` apply the curated
   result.

The edited file already exists, so the cell below runs end-to-end. Re-running the
`suggest` step simply regenerates the suggestions; it does **not** overwrite your
edited file.

> Note: we apply the extensions to the single shared `fault_data` before meshing,
> so the surfaces built in Stage 1 already reflect them.""",
"""### トレース延長（人間によるレビューを含む）

隣接する断層がきちんと接続して綺麗にカットできるよう、断層は走向方向に少し延長
する必要がしばしばあります。これは**2段階・人手レビュー**の処理です：

1. `suggest_trace_extensions(...)` が、*提案された*延長の CSV（および GeoJSON）を
   書き出し、GIS で確認できるようにします。
2. それを編集・承認して `trace_extensions_edited.csv` として保存し、
   `read_trace_extensions(...)` と `apply_trace_extensions()` で精査済みの結果を
   適用します。

編集済みファイルは既に存在するので、下のセルは一通り実行できます。`suggest` の
ステップを再実行しても提案が再生成されるだけで、編集済みファイルを上書きする
ことは**ありません**。

> 注：延長はメッシュ生成の前に、共有する単一の `fault_data` に適用します。
> よってステージ1で生成される面には、すでに延長が反映されています。""")

code(
"""# Step 1: propose extensions (for review in GIS). Safe to re-run.
fault_data.suggest_trace_extensions(
    out_file="suggested_trace_extensions.csv",
    geojson_out_file="suggested_trace_extensions.geojson",
    fit_distance=5.e3, extend_distance=40.e3, proximity_threshold=1.e3)

# Step 2: apply the curated, hand-edited extensions.
fault_data.read_trace_extensions(str(TRACE_EXT_EDITED))
fault_data.apply_trace_extensions()""",
"""# ステップ1：延長を提案する（GIS で確認するため）。再実行しても安全。
fault_data.suggest_trace_extensions(
    out_file="suggested_trace_extensions.csv",
    geojson_out_file="suggested_trace_extensions.geojson",
    fit_distance=5.e3, extend_distance=40.e3, proximity_threshold=1.e3)

# ステップ2：精査済み・手編集の延長を適用する。
fault_data.read_trace_extensions(str(TRACE_EXT_EDITED))
fault_data.apply_trace_extensions()""")

# ---------------------------------------------------------------- STAGE 1
md(
"""## Stage 1 — Build fault surfaces

For each curated fault we draw a stack of **depth contours** (lines of equal
depth, every 500 m down to 32 km), then triangulate the surface spanning those
contours. The result is a raw 3-D fault surface.

If the proper surface-meshing step fails for a fault (it can, for awkward
geometries), we fall back to a simpler contour-based mesh and save the offending
contours to a `failed_contours_*.geojson` for inspection — so one bad fault never
stops the whole run.

Each surface is written three ways: **OBJ** (input to Stage 2), **VTK** (for
quick viewing), and the **contours as GeoJSON**.""",
"""## ステージ1 — 断層面の生成

各精査済み断層について、**深度コンター**（等深線。32 km まで 500 m ごと）を
重ねて描き、それらのコンターにまたがる面を三角形分割します。結果として、生の
3次元断層面が得られます。

本来の面メッシュ生成が（扱いにくい形状などで）失敗した場合は、より単純なコンター
ベースのメッシュにフォールバックし、問題のコンターを `failed_contours_*.geojson`
として保存して確認できるようにします。これにより、1つの不良断層が実行全体を止め
ることはありません。

各面は3つの形式で書き出します：**OBJ**（ステージ2の入力）、**VTK**（簡易表示用）、
および **コンターの GeoJSON**。""")

code(
"""for fault in fault_data.curated_faults:
    print(fault.name)
    try:
        fault.generate_depth_contours(DEPTH_CONTOUR_LEVELS, smoothing=False)
        mesh = fault.mesh_fault_surface(check_mesh=False, resolution=MESH_RESOLUTION)
    except Exception as e:
        # Fallback: simple contour mesh, and save the contours that tripped us up.
        print(f"  Failed to mesh {fault.name}: {e}")
        fault.contours.to_file(f"failed_contours_{fault.name}.geojson", driver="GeoJSON")
        mesh = fault.mesh_simple_contours(DEPTH_CONTOUR_LEVELS)

    fault.contours.to_file(str(CONTOUR_DIR / f"{fault.name}_depth_contours.geojson"), driver="GeoJSON")
    mesh.write(str(VTK_DIR / f"{fault.name}_depth_contours.vtk"))
    mesh.write(str(OBJ_DIR / f"{fault.name}_depth_contours.obj"))""",
"""for fault in fault_data.curated_faults:
    print(fault.name)
    try:
        fault.generate_depth_contours(DEPTH_CONTOUR_LEVELS, smoothing=False)
        mesh = fault.mesh_fault_surface(check_mesh=False, resolution=MESH_RESOLUTION)
    except Exception as e:
        # フォールバック：単純なコンターメッシュ。問題となったコンターを保存。
        print(f"  Failed to mesh {fault.name}: {e}")
        fault.contours.to_file(f"failed_contours_{fault.name}.geojson", driver="GeoJSON")
        mesh = fault.mesh_simple_contours(DEPTH_CONTOUR_LEVELS)

    fault.contours.to_file(str(CONTOUR_DIR / f"{fault.name}_depth_contours.geojson"), driver="GeoJSON")
    mesh.write(str(VTK_DIR / f"{fault.name}_depth_contours.vtk"))
    mesh.write(str(OBJ_DIR / f"{fault.name}_depth_contours.obj"))""")

# ---------------------------------------------------------------- STAGE 2 prep
md(
"""## Stage 2 — Cut the surfaces

Now we truncate faults against each other and against the base depth surface, so
the final meshes don't overlap or extend below the model.

### Stage 2 setup

We load:

- the **depth raster** as a PyVista surface (the lower bound everything is
  trimmed to), and
- the **additional cuts** CSV (extra truncation relationships not captured by the
  hierarchy alone).

There is also an optional inverse, `read_excluded_cuts(...)`, for pairs that must
**never** be cut. The cutting step consults both override lists via
`fault_data.should_cut(...)`: `additional_cuts` is loaded here, while
`excluded_cuts` is left as a commented example (no such pairs in this system).
Both default to empty sets, so to use exclusions elsewhere just uncomment the
`EXCLUDED_CUTS` path and its loader.

Then we read the Stage 1 OBJ surfaces back from disk into each fault's `.mesh`.
Reading from disk (rather than reusing in-memory objects) is deliberate: it keeps
this stage self-contained, so you can re-cut without re-running Stage 1.""",
"""## ステージ2 — 面のカット

次に、断層同士、および基底の深度面に対して断層を切断し、最終メッシュが重なったり
モデルより下に伸びたりしないようにします。

### ステージ2の準備

ここで読み込むもの：

- **深度ラスタ**を PyVista 面として（すべてがこれに合わせてトリミングされる下限）、
- **追加カット** CSV（階層だけでは捉えきれない追加の切断関係）。

逆向きの任意機能として、**絶対にカットしない**ペアを指定する
`read_excluded_cuts(...)` もあります。カット処理では両方の上書きリストを
`fault_data.should_cut(...)` で参照します：`additional_cuts` はここで読み込み、
`excluded_cuts` はコメント例として残しています（この断層システムには該当ペアが
ありません）。どちらも既定では空集合なので、他の用途で除外を使うには
`EXCLUDED_CUTS` のパスとその読み込み行のコメントを外すだけです。

その後、ステージ1の OBJ 面をディスクから各断層の `.mesh` に読み戻します。
（メモリ上のオブジェクトを使い回すのではなく）ディスクから読むのは意図的で、
このステージを自己完結させ、ステージ1を再実行せずに再カットできるようにする
ためです。""")

code(
"""depth_pyvista = read_raster(str(DEPTH_RASTER), use_z=True, out_crs=f"EPSG:{EPSG}")
fault_data.read_additional_cuts(str(ADDITIONAL_CUTS))
# OPTIONAL override list of pairs that must NEVER be cut. Same 2-column,
# header-less CSV format as additional_cuts: each row is `fault_to_cut,cutting_fault`.
# Not needed for this fault system, but uncomment to use it elsewhere:
# fault_data.read_excluded_cuts(str(EXCLUDED_CUTS))

# Load the Stage 1 surfaces back onto each fault.
for fault in fault_data.curated_faults:
    obj_file = OBJ_DIR / f"{fault.name}_depth_contours.obj"
    if obj_file.exists():
        fault.mesh = FaultMesh.from_file(str(obj_file))

# Lookup of faults that actually have a mesh, used during cutting.
cutting_dict = {f.name: f for f in fault_data.curated_faults if f.mesh is not None}
print(f"{len(cutting_dict)} faults have meshes available for cutting")""",
"""depth_pyvista = read_raster(str(DEPTH_RASTER), use_z=True, out_crs=f"EPSG:{EPSG}")
fault_data.read_additional_cuts(str(ADDITIONAL_CUTS))
# 任意：絶対にカットしないペアの上書きリスト。additional_cuts と同じ2列・ヘッダ無し
# の CSV 形式で、各行は `fault_to_cut,cutting_fault`。
# この断層システムでは不要ですが、他で使う場合はコメントを外してください：
# fault_data.read_excluded_cuts(str(EXCLUDED_CUTS))

# ステージ1の面を各断層に読み戻す。
for fault in fault_data.curated_faults:
    obj_file = OBJ_DIR / f"{fault.name}_depth_contours.obj"
    if obj_file.exists():
        fault.mesh = FaultMesh.from_file(str(obj_file))

# 実際にメッシュを持つ断層の参照辞書。カット時に使用。
cutting_dict = {f.name: f for f in fault_data.curated_faults if f.mesh is not None}
print(f"{len(cutting_dict)} faults have meshes available for cutting")""")

# ---------------------------------------------------------------- STAGE 2 cut
md(
"""### Apply the cutting hierarchy

We walk the faults in hierarchy order. For each fault we consider every
**higher-priority** fault. The manual override lists are checked first via
`should_cut(...)` (a pair in `additional_cuts` is always cut, one in
`excluded_cuts` never); otherwise we fall back to `decide_whether_to_cut(...)`,
which asks whether they genuinely intersect and are far enough apart to warrant a
cut (the `higher_meshes` give the context of faults that already truncate the
cutter). If a cut is called for, we build a *cutting fragment* from the cutter and
slice this fault with it, keeping the side nearest the fault's own trace.

Finally each fault is trimmed against the depth surface with `cut_mesh_pv(...)`
and the result written as `<name>_cut.obj` — the input to Stage 3.

`fancy_cutting=True` uses the more robust intersection-following cut.""",
"""### カット階層の適用

断層を階層順にたどります。各断層について、すべての**上位**断層を検討します。
まず手動の上書きリストを `should_cut(...)` で確認し（`additional_cuts` のペアは
必ずカット、`excluded_cuts` のペアは決してカットしない）、該当しなければ
`decide_whether_to_cut(...)` にフォールバックして、実際に交差しておりカットに
値するだけ離れているかを判定します（`higher_meshes` は、カッター側を既に切断
している断層の文脈を与えます）。カットすべき場合は、カッターから*カット用
フラグメント*を作り、この断層を切断して、断層自身のトレースに近い側を残します。

最後に各断層を `cut_mesh_pv(...)` で深度面に対してトリミングし、結果を
`<name>_cut.obj` として書き出します。これがステージ3の入力です。

`fancy_cutting=True` は、より堅牢な交差追従カットを使用します。""")

code(
"""threshold = 0.7  # currently unused by the library; kept for API compatibility

for fault_name in fault_data.cutting_hierarchy:
    print(f"Processing {fault_name} ...")
    fault = fault_data.name_dict[fault_name]
    if fault.mesh is None:
        continue

    # Faults ranked above this one are candidates to cut it.
    faults_to_cut = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(fault_name)]
    for cut_name in faults_to_cut:
        if cut_name not in cutting_dict:
            continue
        higher = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(cut_name)]
        higher_meshes = [cutting_dict[n].mesh for n in higher]
        cutting_mesh = cutting_dict[cut_name].mesh
        # Manual overrides win first. should_cut() matches on the curated fault
        # names (reliable) and returns True/False to force/skip a cut, or None to
        # fall through to the geometric test below. (We match on names here rather
        # than via decide_whether_to_cut's excluded_cuts/additional_cuts kwargs,
        # because mesh names carry a "_depth_contours"/"_cut_by_" suffix that would
        # not match the bare names in the CSVs.) additional_cuts/excluded_cuts
        # default to empty sets, so this is a no-op unless you loaded the CSVs.
        decision = fault_data.should_cut(fault_name, cut_name)
        if decision is None:
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
    print(f"  Done: {fault.name}")""",
"""threshold = 0.7  # 現状ライブラリ側では未使用。API互換のため残置

for fault_name in fault_data.cutting_hierarchy:
    print(f"Processing {fault_name} ...")
    fault = fault_data.name_dict[fault_name]
    if fault.mesh is None:
        continue

    # この断層より上位の断層が、これをカットする候補。
    faults_to_cut = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(fault_name)]
    for cut_name in faults_to_cut:
        if cut_name not in cutting_dict:
            continue
        higher = fault_data.cutting_hierarchy[:fault_data.cutting_hierarchy.index(cut_name)]
        higher_meshes = [cutting_dict[n].mesh for n in higher]
        cutting_mesh = cutting_dict[cut_name].mesh
        # 手動の上書きが最優先です。should_cut() は精査済みの断層名で照合し（確実）、
        # カットを強制/抑止する場合は True/False を、そうでなければ下の幾何判定に
        # 委ねる場合は None を返します。（ここで名前照合するのは、メッシュ名に
        # "_depth_contours"/"_cut_by_" の接尾辞が付き、CSV の素の名前と一致しない
        # ためで、decide_whether_to_cut の excluded_cuts/additional_cuts 引数は
        # 使いません。）両リストは既定で空集合なので、CSV を読み込まない限り
        # 何も起きません。
        decision = fault_data.should_cut(fault_name, cut_name)
        if decision is None:
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

    # 基底の深度面に対してトリミングし、最終カットメッシュを保存。
    try:
        depth_trimmed = fault_data.name_dict[fault.name].mesh.cut_mesh_pv(
            depth_pyvista, fault_trace=fault.original_nztm_trace_array)
        fault_data.name_dict[fault.name].mesh = depth_trimmed
    except Exception as e:
        print(f"  Error trimming depth for {fault.name}: {e}")
        continue
    fault_data.name_dict[fault.name].mesh.mesh.write(str(FINAL_DIR / f"{fault.name}_cut.obj"))
    print(f"  Done: {fault.name}")""")

# ---------------------------------------------------------------- STAGE 3
md(
"""## Stage 3 — Remesh with MMG

The cut surfaces have uneven, often sliver-shaped triangles left over from the
contouring and cutting. This stage cleans each one into **uniform, near-
equilateral triangles** of `TARGET_SIZE` using **MMG** (`mmgs`).

### Why MMG and not a parametrise-and-regenerate remesher?

A reparametrising remesher (e.g. gmsh `classifySurfaces` + `createGeometry`)
flattens the whole surface onto a 2-D plane and re-triangulates there. For a
long, **warped** fault sheet like the Alpine Fault, that flattening *folds over
itself*, and the folded region is silently dropped — **about half of the Alpine
Fault was being deleted.**

MMG instead remeshes the surface **in place**, using local edge
split/collapse/swap and moving nodes back onto the original surface. There is no
global flattening to fold, so:

- nothing is cut off — the fault outline and full extent are preserved
  (area preserved to >99.8 % across all faults), and
- triangle quality jumps from a median of ~0.74 to ~0.97 (1.0 = equilateral).

`HAUSD` controls how tightly the new mesh hugs the original surface; `RIDGE_ANGLE_DEG`
tells MMG which sharp kinks to preserve rather than smooth.""",
"""## ステージ3 — MMG でリメッシュ

カット済みの面には、コンター化やカットの結果として、不揃いでしばしば細長い
（スライバー状の）三角形が残っています。このステージでは、**MMG**（`mmgs`）を
使って各面を `TARGET_SIZE` の**均一・ほぼ正三角形**に整えます。

### なぜ「再パラメータ化して再生成する」リメッシャーではなく MMG なのか？

再パラメータ化型のリメッシャー（例：gmsh の `classifySurfaces` + `createGeometry`）
は、面全体を2次元平面に平坦化し、そこで再三角形分割します。アルパイン断層のような
長く**反った**断層シートでは、この平坦化が*自分自身に折り重なり*、折り重なった
領域が黙って捨てられます。**アルパイン断層の約半分が削除されていました。**

一方 MMG は、局所的な辺の分割・統合・スワップを用い、ノードを元の面上に戻しながら
面を**その場で**リメッシュします。折り重なる原因となる大域的な平坦化が無いため：

- 何も切り落とされず、断層の輪郭と全範囲が保持されます
  （全断層で面積を 99.8 % 超保持）、
- 三角形の品質が中央値 約0.74 から 約0.97 に向上します（1.0 が正三角形）。

`HAUSD` は新しいメッシュが元面にどれだけ密着するかを制御し、`RIDGE_ANGLE_DEG` は
平滑化せずに保持すべき鋭い折れ目を MMG に指示します。""")

code(
"""remeshed_paths = []
for obj_path in sorted(FINAL_DIR.glob("*_cut.obj")):
    name = obj_path.stem.removesuffix("_cut")
    out_path = REMESH_DIR / f"{name}_remeshed.vtk"
    print(f"Remeshing {name}")
    try:
        fault_mesh = FaultMesh.from_file(str(obj_path))
        remeshed = fault_mesh.remesh(target_size=TARGET_SIZE, hausd=HAUSD,
                                     ridge_angle_deg=RIDGE_ANGLE_DEG, verbose=False)
        meshio.write(str(out_path), remeshed.mesh)
        remeshed_paths.append(out_path)
    except Exception as e:
        print(f"  failed: {e}")""",
"""remeshed_paths = []
for obj_path in sorted(FINAL_DIR.glob("*_cut.obj")):
    name = obj_path.stem.removesuffix("_cut")
    out_path = REMESH_DIR / f"{name}_remeshed.vtk"
    print(f"Remeshing {name}")
    try:
        fault_mesh = FaultMesh.from_file(str(obj_path))
        remeshed = fault_mesh.remesh(target_size=TARGET_SIZE, hausd=HAUSD,
                                     ridge_angle_deg=RIDGE_ANGLE_DEG, verbose=False)
        meshio.write(str(out_path), remeshed.mesh)
        remeshed_paths.append(out_path)
    except Exception as e:
        print(f"  failed: {e}")""")

md(
"""### Merge for inspection

Finally we merge every remeshed surface into one VTK so the whole fault network
can be opened and checked in one go (e.g. in ParaView).""",
"""### 確認用にマージ

最後に、リメッシュ済みの全ての面を1つの VTK にマージし、断層ネットワーク全体を
一度に開いて確認できるようにします（例：ParaView）。""")

code(
"""meshes = [pv.from_meshio(meshio.read(str(p))) for p in remeshed_paths]
combined = pv.merge(meshes)
combined.save(str(MERGED_VTK))
print(f"Merged {len(meshes)} remeshed surfaces -> {MERGED_VTK}")""",
"""meshes = [pv.from_meshio(meshio.read(str(p))) for p in remeshed_paths]
combined = pv.merge(meshes)
combined.save(str(MERGED_VTK))
print(f"Merged {len(meshes)} remeshed surfaces -> {MERGED_VTK}")""")

# ---------------------------------------------------------------- viz
md(
"""### (Optional) Quick 3-D view

A quick look at the merged result. This needs an interactive/display-capable
environment; if it doesn't render inline, set a PyVista Jupyter backend
(e.g. `pv.set_jupyter_backend("trame")`) or open `merged_mesh_remeshed.vtk` in
ParaView instead.""",
"""### （任意）簡易3次元表示

マージ結果をざっと確認します。これは対話的・表示可能な環境が必要です。インライン
表示されない場合は、PyVista の Jupyter バックエンドを設定するか
（例：`pv.set_jupyter_backend("trame")`）、`merged_mesh_remeshed.vtk` を
ParaView で開いてください。""")

code(
"""combined = pv.read(str(MERGED_VTK))
plotter = pv.Plotter()
plotter.add_mesh(combined, show_edges=True, color="lightgray")
plotter.show()""",
"""combined = pv.read(str(MERGED_VTK))
plotter = pv.Plotter()
plotter.add_mesh(combined, show_edges=True, color="lightgray")
plotter.show()""")


def build(lang_index, path):
    nb = new_notebook()
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3 (leapfrog-fault-models)",
                       "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    for kind, en, ja in CELLS:
        text = (en, ja)[lang_index]
        nb.cells.append(new_markdown_cell(text) if kind == "md" else new_code_cell(text))
    with open(path, "w", encoding="utf-8") as f:
        nbf.write(nb, f)
    print("wrote", path)


build(0, "central_nz_pipeline.ipynb")
build(1, "central_nz_pipeline_ja.ipynb")
