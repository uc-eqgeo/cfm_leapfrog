#!/bin/bash
bounds="-R412000/2828000/3944000/6500000"


Df_file="ESupp1_RuptureDepthSurfaces_v1.0.csv"

#awk 'BEGIN {FS=","} {if(NR>1) print $1,$2,$16}' ${Df_file} |gmt blockmean -I2000/2000 ${bounds} > df_means.xyz
#gmt surface df_means.xyz -GDf.grd ${bounds} -I2000/2000
#gmt grdclip Df.grd -SrNaN/0 -GDf_no_nans.grd
#gmt grdmath Df_no_nans.grd -1000 MUL = Df_no_smoothing_km.grd
#gmt grdfilter Df_no_smoothing_km.grd -D0 -Fb50000 -GDf_smoothed_50k.grd
#gdal_translate -a_srs epsg:2193 Df_smoothed_50k.grd Df_smoothed_50k.tif

### Distance over which smoothing should occur
buffer_size=2.5e4
### Distance to shift subduction zone surfaces by
shift_size=500.

#
#### Hikurangi ##########################################################################################################
#### Fit surface through Hikurangi vertices
#gmt surface hikurangi_points.csv -Ghikurangi.grd ${bounds} -I2000/2000
#### Shift Hikurangi surface upwards by 500 m
#gmt grdmath hikurangi.grd ${shift_size} ADD = hikurangi_shifted.grd
#### Create mask of Hikurangi surface using outline csv (inside set to 1, outside set to 0)
#gmt grdmask hikurangi_outline.csv -Ghikurangi_mask_inside.grd ${bounds} -I2000/2000 -N0/1/1
#### Create mask of Hikurangi surface using outline csv (inside set to 0, outside set to 1)
#gmt grdmath hikurangi_mask_inside.grd NOT = hikurangi_mask_outside.grd
#### Create buffered grid of Hikurangi outline using buffered outline csv (inside set to 1, outside set to 0)
#gmt grdmask hikurangi_outline_buffered.csv -Ghikurangi_mask_inside_buffer.grd ${bounds} -I2000/2000 -N0/1/1
#### Create mask of Hikurangi surface using buffered outline csv (inside set to 0, outside set to 1)
#gmt grdmath hikurangi_mask_inside_buffer.grd NOT = hikurangi_mask_outside_buffer.grd
#### Trim maximum surface to be within buffer (not that important first time, but could be for e.g. Puysegur)
#gmt grdmath hikurangi_shifted.grd hikurangi_mask_inside_buffer.grd MUL = hikurangi_shifted_inside_buffer.grd
#### Set depths outside buffer so they are definitely deeper than Dfc surface (for later comparisons)
#gmt grdmath hikurangi_mask_outside_buffer.grd -80e3 MUL = hikurangi_depths_outside_buffer.grd
#### Combine depths inside and outside buffer
#gmt grdmath hikurangi_shifted_inside_buffer.grd hikurangi_depths_outside_buffer.grd ADD = hikurangi_modified_depths.grd
#### Find which is shallower at each point, the Hikurangi surface or the Df surface
#gmt grdmath hikurangi_modified_depths.grd Df_smoothed_50k.grd MAX = hikurangi_max.grd
#### Find distance of each point from outline CSV (and normalise to multiple of buffer distance)
#gmt grdmath hikurangi_max.grd hikurangi_outline.csv LDIST ${buffer_size} DIV = hik_distance.grd
#### Set all distances > 1 to 1.
#gmt grdclip hik_distance.grd -Sa1/1 -Ghik_distance_normalized.grd
#### Multiply distance by Hikurangi surface mask (so that distance is 0 at Hikurangi surface and 1 at buffer edge)
#gmt grdmath hik_distance_normalized.grd hikurangi_mask_outside.grd MUL = hikurangi_for_blend.grd
#### Linear blending of Hikurangi max surface and Df surface
#gmt grdmath Df_smoothed_50k.grd hikurangi_max.grd SUB hikurangi_for_blend.grd MUL hikurangi_max.grd ADD = hikurangi_blended.grd
#
#### Puysegur ###########################################################################################################
## Create surface for Puysegur from vertices
#gmt surface puysegur_points.txt -Gpuysegur.grd ${bounds} -I2000/2000
## Shift Puysegur surface upwards by 500 m
#gmt grdmath puysegur.grd ${shift_size} ADD = puysegur_shifted.grd
## Create mask of Puysegur surface using outline csv (inside set to 1, outside set to 0)
#gmt grdmask puysegur_outline.gmt -Gpuysegur_mask_inside.grd ${bounds} -I2000/2000 -N0/1/1
## Create mask of Puysegur surface using outline csv (inside set to 0, outside set to 1)
#gmt grdmath puysegur_mask_inside.grd NOT = puysegur_mask_outside.grd
## Create buffered grid of Puysegur outline using buffered outline csv (inside set to 1, outside set to 0)
#gmt grdmask puysegur_outline_buffered.gmt -Gpuysegur_mask_inside_buffer.grd ${bounds} -I2000/2000 -N0/1/1
## Create mask of Puysegur surface using buffered outline csv (inside set to 0, outside set to 1)
#gmt grdmath puysegur_mask_inside_buffer.grd NOT = puysegur_mask_outside_buffer.grd
## Trim maximum surface to be within buffer (not that important first time, but could be for e.g. Puysegur)
#gmt grdmath puysegur_shifted.grd puysegur_mask_inside_buffer.grd MUL = puysegur_shifted_inside_buffer.grd
## Set depths outside buffer so they are definitely deeper than Dfc surface (for later comparisons)
#gmt grdmath puysegur_mask_outside_buffer.grd -80e3 MUL = puysegur_depths_outside_buffer.grd
## Combine depths inside and outside buffer
#gmt grdmath puysegur_shifted_inside_buffer.grd puysegur_depths_outside_buffer.grd ADD = puysegur_modified_depths.grd
## Find which is shallower at each point, the Puysegur surface or the Df surface with Hikurangi blended in
#gmt grdmath puysegur_modified_depths.grd hikurangi_blended.grd MAX = puysegur_max.grd
## Find distance of each point from outline CSV (and normalise to multiple of buffer distance)
#gmt grdmath puysegur_max.grd puysegur_outline.gmt LDIST ${buffer_size} DIV = puy_distance.grd
## Set all distances > 1 to 1.
#gmt grdclip puy_distance.grd -Sa1/1 -Gpuy_distance_normalized.grd
## Multiply distance by Puysegur surface mask (so that distance is 0 at Puysegur surface and 1 at buffer edge)
#gmt grdmath puy_distance_normalized.grd puysegur_mask_outside.grd MUL = puysegur_for_blend.grd
## Linear blending of Puysegur max surface and Df surface with Hikurangi blended in
#gmt grdmath hikurangi_blended.grd puysegur_max.grd SUB puysegur_for_blend.grd MUL puysegur_max.grd ADD = puysegur_blended.grd
#
#gdal_translate -a_srs epsg:2193 puysegur_blended.grd puysegur_blended.tif

### Add in local modifications for Hannu
# Mask out the areas that will be changed
gmt grdmask TectonicDomains_NZ_CFM_v1_0_HikurangiOR_modified.gmt TectonicDomains_NZ_CFM_v1_0_Puysegur_CaswellHighOR_modified.gmt Fiordland_Dfc_extent.gmt -Ghannu_mods_outside.grd ${bounds} -I2000/2000 -N1/0/0
# Create mask to isolate areas of Hikurangi OR that will be changed and multiply by -15000
gmt grdmask TectonicDomains_NZ_CFM_v1_0_HikurangiOR_modified.gmt -Ghikurangi_or_mask_inside.grd ${bounds} -I2000/2000 -N0/1/1
gmt grdmath hikurangi_or_mask_inside.grd -15000.0 MUL = hikurangi_or_replacement.grd
# Create mask to isolate areas of Puysegur OR that will be changed and multiply by -6000
gmt grdmask TectonicDomains_NZ_CFM_v1_0_Puysegur_CaswellHighOR_modified.gmt -Gpuysegur_or_mask_inside.grd ${bounds} -I2000/2000 -N0/1/1
gmt grdmath puysegur_or_mask_inside.grd -6000.0 MUL = puysegur_or_replacement.grd
# Create mask to isolate areas of Fiordland Dfc that will be changed and multiply by Dfc surface
gmt grdmask Fiordland_Dfc_extent.gmt -Gfiordland_dfc_mask_inside.grd ${bounds} -I2000/2000 -N0/1/1
gmt grdmath fiordland_dfc_mask_inside.grd Df_smoothed_50k.grd MUL = fiordland_dfc_replacement.grd

# Set areas that will be changed to 0 in unedited grid
gmt grdmath puysegur_blended.grd hannu_mods_outside.grd MUL = puysegur_for_mods.grd
# Add in the modified areas
gmt grdmath puysegur_for_mods.grd puysegur_or_replacement.grd ADD hikurangi_or_replacement.grd ADD fiordland_dfc_replacement.grd ADD = with_hannu_mods.grd

gdal_translate -a_srs epsg:2193 with_hannu_mods.grd with_hannu_mods.tif

