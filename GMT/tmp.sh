#!/usr/bin/env bash


####export GMT_TEXT_FONT=SimSun

gmt begin FOCAL png

gmt set FONT 6p,Times-Roman
gmt set MAP_FRAME_PEN 0.4p,black
gmt set MAP_FRAME_AXES  WNes
gmt set MAP_FRAME_TYPE = plain 
gmt set FORMAT_GEO_MAP=ddd:mm:ssF
##gmt set  MAP_ANNOT_MIN_SPACING = plain 15





# Base plot
gmt grdimage @earth_relief_30s -JM9c -t30 -Cworld -I+d   -R120.5/122.25/22.5/24.5 -Ba0.5f0.25
gmt colorbar -Bxaf+l"Elevation (m)" -DJML+w4c+o1.2c/0c+ml -Bxa1000f 



# National boundary
gmt plot CN-faults.gmt -W0.77p,gray58
####gmt plot mfa.gmt -W0.77p,gray58
#gmt plot CN-block-L1.gmt -W1.1p,gray
#gmt plot CN-block-L2.gmt -W1.1p,gray



cat << EOF > Icpt.cpt
 0   0-1-1   20   0-1-1
20  60-1-1   40  60-1-1
40 120-1-1   60 120-1-1
60 240-1-1  100 240-1-1
EOF

gmt makecpt -Cwhite,yellow,red   -T0/40/1  -H -D -Z > Icpt.cpt
gmt colorbar -CIcpt.cpt -DjBR+w1.25c/0.15c+o0.8c/0.25c+ml -Bx+l'source depth' -By+lkm -G0/40
gmt meca tmp.focal  -Ba -A+s0.2c -Sa0.35c -CIcpt.cpt 










gmt end show






