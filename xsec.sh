#!/bin/bash
xmin=`sort -k 1 -n ./TongaAttenData/Model_xsec1.dat | head -1 | awk -F ' ' '{print $1}'`
xmax=`sort -k 1 -n ./TongaAttenData/Model_xsec1.dat | tail -1 | awk -F ' ' '{print $1}'`
ymin=`sort -k 2 -n ./TongaAttenData/Model_xsec1.dat | head -1 | awk -F ' ' '{print $2}'`
ymax=`sort -k 2 -n ./TongaAttenData/Model_xsec1.dat | tail -1 | awk -F ' ' '{print $2}'`
R=0/$xmax/$ymin/$ymax
J="X-20c/-14c" 


# end2=( 178.10/-16.32 177.60/-17.32 177.10/-18.32 )
# end1=( -173.90/-20.31 -174.40/-21.31 -174.90/-22.31 )

end2=( 178.10/-16.32 177.60/-17.32 177.10/-18.32 )
end1=( 186.1/-20.31 185.6/-21.31 185.1/-22.31 )

abs() {
  local number=$1
  echo $((number < 0 ? -number : number))
}

# Example usage:
number=-5
#absolute_value=$(abs $number)
#echo "The absolute value of $number is $absolute_value"

# gmt begin xsec pdf A+m0.05c
# gmt basemap -R$R -J$J -Bxaf -Byaf+l"Depth (km)" -BWSen 
# gmt makecpt -CattenP2.cpt -T0/20/0.1 
# awk -F ' ' '{print $1, $2, $3}' input.dat | gmt plot -Ss0.5c -t10 -C

for xsec in 1 2 3
do
    gmt begin xsec$xsec pdf A+m0.05c
    gmt basemap -R$R -J$J -Bxa100f -Bya100f+l"Depth (km)" -BWSen 
    gmt makecpt -Cjet -D -T0/20/0.1 

    gmt triangulate ./TongaAttenData/Mask_xsec$xsec.dat -R$R -I1/1 -Goutput.grd
    gmt grdimage output.grd -C

    gmt triangulate ./TongaAttenData/STD_xsec$xsec.dat -I0.01/0.01 -Goutput.contour
    gmt grdcontour output.contour -C1 -A1 -Wthinnest,black 
    
    gmt plot ./TongaAttenData/EQ_xsec$xsec.dat -Sc0.3c -W1p,white
    gmt colorbar -DJMR+w8c+o0.8c/0c -Bxa3+l"1000/Qp"

    #Plot topo on top
    gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} -G0.1 > gc.xy
    gmt grdtrack gc.xy -G@earth_relief_02m.grd  > topo.dat
    # paste gc.xyt slice.gc | awk ' { print 111*$3, r, $6 } ' r=$dep | cat >> slice.out
    awk '{print 111*$3, -1*$4}' topo.dat > xy.dat
    tmin=`sort -k 1 -n xy.dat | head -1 | awk -F ' ' '{print $1}'`
    tmax=`sort -k 1 -n xy.dat | tail -1 | awk -F ' ' '{print $1}'`
    gmt basemap -R$tmin/$tmax/-1000/5000 -JX20c/-1.5c -Bx+af+I -Byaf -Bwsen -Y14c
    gmt plot xy.dat -W3p
    if (($xsec == 1)); then
        echo $tmin -1000 "A" | gmt text -F+f16p -N -D0c/0.5c
        echo $tmax -1000 "A'" | gmt text -F+f16p -N -D0c/0.5c
    elif (($xsec == 2)); then
        echo $tmin -1000 "B" | gmt text -F+f16p -N -D0c/0.5c
        echo $tmax -1000 "B'" | gmt text -F+f16p -N -D0c/0.5c
    elif (($xsec == 3)); then
        echo $tmin -1000 "C" | gmt text -F+f16p -N -D0c/0.5c
        echo $tmax -1000 "C'" | gmt text -F+f16p -N -D0c/0.5c
    fi
    
    if (($xsec == 1)); then
        CLSC_projection=`cat gmtdata/CLSC.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$CLSC_projection" ]; then
            echo "too far away"
        else
            echo $CLSC_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $CLSC_projection -1000 "CLSC" | gmt text -F+f16p -N -D-0.45c/1.0c
        fi

        ELSC_projection=`cat gmtdata/ELSC.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$ELSC_projection" ]; then
            echo "too far away"
        else
            echo $ELSC_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $ELSC_projection -1000 "ELSC" | gmt text -F+f16p -N -D0.45c/1.0c
        fi

        VFR_projection=`cat gmtdata/VFR.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$VFR_projection" ]; then
            echo "too far away"
        else
            echo $VFR_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $VFR_projection -1000 "VFR" | gmt text -F+f16p -N -D0c/1.0c
        fi
        
        Tofua_projection=`awk -F ' ' '{if(-173>=$2 && $2>=-177) print $2, $1}' volcanoes.lst | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        echo $Tofua_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
        echo $Tofua_projection -1000 "Tofua" | gmt text -F+f16p -N -D0c/1.0c
    else
        CLSC_projection=`cat gmtdata/CLSC.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$CLSC_projection" ]; then
            echo "too far away"
        else
            echo $CLSC_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $CLSC_projection -1000 "CLSC" | gmt text -F+f16p -N -D0c/1.0c
        fi

        ELSC_projection=`cat gmtdata/ELSC.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$ELSC_projection" ]; then
            echo "too far away"
        else
            echo $ELSC_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $ELSC_projection -1000 "ELSC" | gmt text -F+f16p -N -D0c/1.0c
        fi

        VFR_projection=`cat gmtdata/VFR.dat | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        if [ -z "$VFR_projection" ]; then
            echo "too far away"
        else
            echo $VFR_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
            echo $VFR_projection -1000 "VFR" | gmt text -F+f16p -N -D0c/1.0c
        fi

        # Tofua_projection=`awk -F ' ' '{print $2, $1}' volcanoes.lst | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(0.5>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        Tofua_projection=`awk -F ' ' '{if(-173>=$2 && $2>=-177) print $2, $1}' volcanoes.lst | gmt project -C${end2[$((xsec-1))]} -E${end1[$((xsec-1))]} | awk -F ' ' '{if(1.0>=$4 && $4>=0) print}' | sort -k 4 -n | head -1 | awk -F ' ' '{print $3*111}'`
        echo $Tofua_projection -3500 90 0.65c | gmt plot -Sv0.5c+e -W2p -Gblack -N 
        echo $Tofua_projection -1000 "Tofua" | gmt text -F+f16p -N -D0.2c/1.0c
        echo $Tofua_projection -1000 "Tofua"

    fi

    # gmt grdimage -R$region -J$projection -Cgray @earth_relief_01m -I+d
    # gmt makecpt -Cjet -T0.5/20/0.1 
    # # awk -F ' ' '{print $1, $2, $3}' ./TongamapNew/TongaMap$depth.txt | gmt plot -Ss0.05c -t10 -C
    # #gmt xyz2grd ./TongaAttenData/Tonga_Map_Mask_$depth.txt -I0.1/0.1 -Goutput.grd
    # gmt psclip triangle.coordinate
    # # gmt triangulate ./TongaAttenData/Tonga_Map_Transparency_$depth.txt -I0.01/0.01 -Goutput.alpha

    # gmt triangulate ./TongaAttenData/Tonga_Map_Mask_$depth.txt -I0.01/0.01 -Goutput.grd
    # # gmt grdimage output.grd -C -Ioutput.alpha
    # gmt grdimage output.grd -C 
    # gmt triangulate ./TongaAttenData/Tonga_Map_Uncertainty_$depth.txt -I0.01/0.01 -Goutput.contour
    # gmt grdcontour output.contour -C1 -A2 -Wthinnest,black 

    # gmt psclip -C
    
    # gmt plot -A -W2.5p ./gmtdata/lau_neovolcanic.xy

    # echo 177 -23 "$depth km" | gmt text -Gwhite -F+f25p -C15%/15%
    # gmt colorbar -DJMR+w8c+o0.8c/0c -Bxa3+l"1000/Qp"
    # gmt coast -R$region -J$projection -Dh0.5p -W1.5p 
    gmt end



rm gc.xy output.grd topo.dat xy.dat output.contour
done
