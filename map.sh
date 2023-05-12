#!/bin/bash

region="175/188/-24/-14"
mapregion="175.78/187.12/-23.73/-14.80"
projection="M16c"
# gmt set FONT_LABEL 25p 
gmt set FONT_ANNOT_PRIMARY 20p
gmt set FONT_ANNOT_SECONDARY 30p
for depth in "50" "100" "300" "500"
do
    gmt begin map$depth pdf A+m0.5c
    gmt basemap -R$region -J$projection -Baf -BWSen
    gmt grdimage -R$region -J$projection -Cgray @earth_relief_01m -I+d
    gmt makecpt -Cjet -D -T0/20/0.1 
    # awk -F ' ' '{print $1, $2, $3}' ./TongamapNew/TongaMap$depth.txt | gmt plot -Ss0.05c -t10 -C
    #gmt xyz2grd ./TongaAttenData/Tonga_Map_Mask_$depth.txt -I0.1/0.1 -Goutput.grd
    gmt psclip rectangular.coordinate
    # gmt triangulate ./TongaAttenData/Tonga_Map_Transparency_$depth.txt -I0.01/0.01 -Goutput.alpha

    gmt triangulate ./TongaAttenData/Tonga_Map_Mask_$depth.txt -I0.01/0.01 -Goutput.grd
    # gmt grdimage output.grd -C -Ioutput.alpha
    gmt grdimage output.grd -C 
    gmt triangulate ./TongaAttenData/Tonga_Map_Uncertainty_$depth.txt -I0.01/0.01 -Goutput.contour
    gmt grdcontour output.contour -C1 -A2 -Wthinnest,black 

    gmt psclip -C
    
    gmt plot -A -W2.5p ./gmtdata/lau_neovolcanic.xy

######Plot the 3 lines#######    
gmt pstext -Gwhite -F+f16,Helvetica-Bold,black+jMC -N << EOF 
178.00 -16.20 A
-174.00 -20.50 A'
177.50 -17.20 B
-174.50 -21.50 B'
177.00 -18.20 C
-175.00 -22.50 C'
EOF

gmt psxy -W2p <<EOF 
178.10      -16.32
-173.90     -20.31
>
177.60      -17.32
-174.40     -21.31
>
177.10      -18.32
-174.90     -22.31
EOF
############################


# gmt psxy -W2p  << EOF
# 176.88 -16.95
# 186.04 -21.53
# EOF


###############################
    echo 177 -23 "$depth km" | gmt text -Gwhite -F+f25p -C15%/15%
    gmt colorbar -DJMR+w8c+o0.8c/0c -Bxa3+l"1000/Qp"
    gmt coast -R$region -J$projection -Dh0.5p -W1.5p 
    gmt end



rm output.grd output.contour
done







# #     outf2 = "map"$depth
        
#     # TITLE=$depth"  km"

#     # gmt begin $outf2 pdf A+m0.5c
# gmt begin map50 pdf A+m0.5c
# gmt basemap -R$region -J$projection -Baf -BWSen
# # +t"Attenuation in map view (50 km)" 
# gmt grdimage -R$region -J$projection -Cgray @earth_relief_01m -I+d
# # gmt coast -R$region -J$projection 

# gmt plot -A -W1p  << EOF
# 175.78 -19.14
# 177.95 -14.80 
# 187.12 -19.38
# 184.95 -23.73
# 175.78 -19.14
# EOF

# # gmt surface TongaMap50.txt -R$mapregion -I0.2/0.2 -Goutput.nc

# # gmt grd2cpt @output.nc -T0.5/15/1 -Cjet 

# # gmt grdimage @output.nc -C
# gmt makecpt -Cjet -T0.5/20/0.1 
# awk -F ' ' '{print $1, $2, $3}' TongaMap50.txt | gmt plot -Ss0.05c -t10 -C
# # gmt plot -A -W2.5p ./gmtdata/CLSC.dat
# # gmt plot -A -W2.5p ./gmtdata/ELSC.dat
# # gmt plot -A -W2.5p ./gmtdata/VFR.dat
# gmt plot -A -W2.5p ./gmtdata/lau_neovolcanic.xy

# gmt coast -R$region -J$projection -Dh0.5p -W0.5p 
# gmt pstext -Gwhite -F+f16,Helvetica-Bold,black+jMC -N << EOF 
# 178.00 -16.20 A
# -174.00 -20.50 A'
# 177.50 -17.20 B
# -174.50 -21.50 B'
# 177.00 -18.20 C
# -175.00 -22.50 C'
# EOF

# gmt psxy  -W1 <<EOF 
# 178.10      -16.32
# -173.90     -20.31
# >
# 177.60      -17.32
# -174.40     -21.31
# >
# 177.10      -18.32
# -174.90     -22.31
# EOF

# echo 177 -23 "50 km" | gmt text -Gwhite -F+f25p -C15%/15%
# gmt colorbar -DJMR+w8c+o0.8c/0c -Bxa3+l"1000/Qp"
# gmt end


# done

