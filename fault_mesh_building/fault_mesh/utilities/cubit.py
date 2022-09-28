"""
Functions to make journal files for remeshing using cubit.
"""
import os
import geopandas as gpd
from shapely.geometry import LineString, Polygon


def make_journal_file_commands(line: LineString, outfile, outmesh: str, top_z=0., depth=2.e4):
    out_str = ""
    for coord in list(line.coords):
        x, y = coord[:2]
        out_str += f"create vertex location {x:.4f} {y:.4f} {top_z:.2f}\n"

    vertex_numbers = ",".join([str(i + 1) for i in range(len(list(line.coords)) - 1)])
    surface_numbers = ",".join([str(i + 1) for i in range(len(list(line.coords)) + 1)])
    out_str += f"create surface vertex {vertex_numbers}\n"
    out_str += f"sweep surface 1 vector 0 0 -1 distance {depth:.2f}\n"
    out_str += f"surface {surface_numbers} scheme trimesh geometry approximation angle 15\n"
    out_str += f"trimesher surface gradation 1.3\n"
    out_str += f"trimesher geometry sizing on\n"

    out_str += f"mesh surface {surface_numbers}\n"
    out_str += f"""export stl "{outmesh}" mesh overwrite\n"""
    out_str += f"""exit()\n"""

    with open(outfile, "w") as out_id:
        out_id.write(out_str)


def make_journal_file_polygon(polygon: Polygon, outjou: str, outmesh: str, top_z=0., depth=2.e4):
    return make_journal_file_commands(polygon.exterior, outjou, outmesh, top_z, depth)


def make_journal_file_multi(gis_file: str, out_directory: str, top_z=0., depth=2.e4):
    assert os.path.exists(gis_file)
    data = gpd.read_file(gis_file).explode()
    assert "fault_name" in data.columns
    assert all([isinstance(x, Polygon) for x in list(data.geometry)])

    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    for i, row in data.iterrows():
        out_name = f"{out_directory}/{row.fault_name}.jou"
        out_mesh_name = f"{row.fault_name}_volume.stl"
        make_journal_file_polygon(row.geometry, out_name, out_mesh_name, top_z, depth)
