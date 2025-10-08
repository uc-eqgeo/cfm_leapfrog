#!/bin/bash
# To create GMT grid from xyz with no smoothing
gmt surface d95.d0595.3d.050.lst.rslv.lt.20km.Mge1.5.W-0.20-0.05-0.05-120-150-20-50.ge20 -I0.05 -T0.25 -Gjapan_no_smoothing.grd -R120/150/20/50
# Smooth using 20 km Gaussian (comment and change subsequent file names if you don't want to smooth)
gmt grdfilter japan_no_smoothing.grd -D3 -Fg20 -Gjapan_smoothed.grd
# Multiply by -1000 to get positive up and metre depths
gmt grdmath japan_smoothed.grd -1000 MUL = japan_smoothed_metres.grd
# Convert to tif
gmt grdconvert japan_smoothed_metres.grd japan_smoothed_metres.tif=gd:GTiff
# convert to UTM: EPSG:32654
gdalwarp -s_srs EPSG:4326 -t_srs EPSG:32654 japan_smoothed_metres.tif japan_smoothed_metres_UTM.tif