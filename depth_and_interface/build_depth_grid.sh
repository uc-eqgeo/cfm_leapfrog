#!/bin/bash
bounds=-R572000/2443000/4223000/6202000
#gmt grd2xyz hikurangi_interface.tif -s > hikurangi_interface.xyz
#gmt surface hikurangi_interface.xyz -Ghikurangi_expanded.grd -M2000 -I2000 $bounds

gmt grdclip hikurangi_expanded.grd -SrNaN/0 -Ghikurangi_expanded_no_nans.grd

gmt grdmask depth_polygons.gmt -Gdepth_polygons.grd -NZ -Rhikurangi_expanded.grd
gmt grdmask depth_polygons.gmt -Gclipping_region.grd -N1/0/0 -Rhikurangi_expanded.grd

gmt grdmath clipping_region.grd hikurangi_expanded_no_nans.grd MUL 1000 MUL depth_polygons.grd ADD = depths_no_smoothing.grd
gmt grdfilter depths_no_smoothing.grd -D0 -Fb30000 -Gdepths_smoothed.grd

gdal_translate -a_srs epsg:2193 depths_smoothed.grd depths_smoothed.tif

