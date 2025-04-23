import numpy as np
import pyvista as pv
import rioxarray
from pyproj import transform

def read_raster(filename, out_crs="EPSG:2193", use_z=False):
    """Read a raster to a ``pyvista.StructuredGrid``.

    This will handle coordinate transformations.
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
    lon, lat = transform(data.rio.crs, out_crs, pts[:, 0], pts[:, 1])
    mesh.points[:, 0] = lon
    mesh.points[:, 1] = lat
    mesh["data"] = values.reshape(mesh.n_points, -1, order="F")
    return mesh