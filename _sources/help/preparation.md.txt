# Preparing inputs for use

## Shapefile fields
The primary input to the workflow is a shapefile of similar (GeoJSON and GPKG also work).
This shapefile includes the linework for surface traces, but should also include enough metadata to unambiguously constrain the fault geometry.

A minimal example of a shapefile with the required fields is provided in `tutorial_gis/central_nz_minimal_data.shp`.

The required fields are:
- `Name`: a unique identifier for each fault
- `Fault_ID`: a unique identifier (number) for each fault
- `Dip_pref`: the preferred dip of the fault in degrees
- `Dip_dir`: the dip direction of the fault in compass direction (e.g. "N", "NE", "E", "SE", "S", "SW", "W", "NW"). Note that two possible dip azimuths will be calculated based on the fault trace. Inclusion of this letter-based field allows the workflow to choose between these two possible azimuths, and also allows consistency checks (e.g., an error will be thrown if "E" is provided for an EW-striking fault).
- `SR_pref`: the preferred slip rate of the fault in mm/yr

The `SR_pref` field is optional, but if it is not provided, the automated method to estimate which faults terminate against each other (cutting hierarchy) will not work.

**Note**: Other fields can be included in the shapefile, but they will be ignored by the workflow.

## Vertex spacing

The workflow requires that the fault traces are represented by a series of vertices. The spacing of these vertices is important for the accuracy of the calculations. We suggest a spacing of 1 km, but this can be adjusted based on the size of fault networks and complexity of the fault geometries. It would be possible to set this point spacing programmatically using python, but in our experience this iterative process is best carried out interatively in GIS software, for example using the [Densify by interval](https://docs.qgis.org/3.34/en/docs/user_manual/processing_algs/qgis/vectorgeometry.html#densify-by-interval) tool in QGIS.  

