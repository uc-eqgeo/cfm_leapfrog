import numpy as np
import pyvista as pv
import rioxarray
from pyproj import Transformer


def read_raster(filename, out_crs="EPSG:2193", use_z=False):
    """Read a raster to a ``pyvista.Surface``.

    This will handle coordinate transformations.

    Adapted from code by Bane Sullivan (PyVista). See
    https://github.com/banesullivan/banesullivan/blob/5e88f9e406cdf12c7c9652071687919d241f7bf2/pyvista-examples/raster.py#L34-L41
    """
    # Read in the data
    data = rioxarray.open_rasterio(filename)
    values = np.asarray(data)
    data.rio.nodata
    nans = values == data.rio.nodata
    if np.any(nans):
        # values = np.ma.masked_where(nans, values)
        values[nans] = np.nan
    # Make a mesh
    xx, yy = np.meshgrid(data["x"], data["y"])
    if use_z and values.shape[0] == 1:
        # will make z-comp the values in the file
        zz = values.reshape(xx.shape)
    else:
        # or this will make it flat
        zz = np.zeros_like(xx)
    mesh = pv.StructuredGrid(xx, yy, zz)
    pts = mesh.points
    transformer = Transformer.from_crs(data.rio.crs, out_crs, always_xy=True)
    lon, lat = transformer.transform(pts[:, 0], pts[:, 1])
    mesh.points[:, 0] = lon
    mesh.points[:, 1] = lat
    mesh["data"] = values.reshape(mesh.n_points, -1, order="F")
    return mesh.extract_surface(algorithm="dataset_surface")