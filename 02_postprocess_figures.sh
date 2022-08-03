#!/bin/bash


echo "Making output directory 'final_figures'..."
mkdir final_figures


echo "Making final figures... "


# Merging two sub-figures
#convert \( figures/figure-stommel_bifurcation_phase.png \) \
#    \( figures/figure-ZATOM_bifurcation_phase.png \) -gravity center -append \
#     figures/merged-phase-diagram.png

#convert \( figures/figure-stommel_bifurcation_phase.png \) \
#    \( figures/figure-ZATOM_bifurcation_phase.png \) -gravity center -append \
#     figures/merged-phase-diagram.png

echo "Doing merging : analytical extended two-box model"
convert \( figures/figure-stommel_bifurcation_analytical_p.png \) \
    \( figures/figure-stommel_bifurcation_analytical_xi.png \) -gravity North +append \
     figures/merged-stommel_bifurcation_analytical.png

#echo "Doing merging : regime diagrams"
#convert \( figures/figure-ZATOM_bifurcation_phase.png \) \
#    \( figures/figure-stommel_bifurcation_phase.png \) -gravity North +append \
#     figures/merged-regime_diagrams.png

echo "Doing merging : freshwater forcing and cartoon"
convert \( figures/figure-forcing.png \) \
    \( figures/cartoon_ZATOM_forcing.png -resize 80% \) -gravity center +append \
     figures/merged-forcing.png


cp manual_figures/* figures/

name_pairs=(
    ZATOM_design.png                              fig01.png 
    merged-forcing.png                            fig02.png
    MOC_dynaimcs.png                              fig03.png
    ZATOM_bifur_analysis_xi.png                   fig04.png
    ZATOM_bifur_analysis_xi_marks.png             fig05.png
    ZATOM_bifur_gamma_xi.png                      fig06.png
    figure-reduced_stommel.png                    fig07.png
    merged-stommel_bifurcation_analytical.png     fig08.png
    figure-stommel_bifurcation_phase.png          fig09.png
    regime_diagrams_comparison.png                fig10.png
    cartoon_dijwei2013_forcing.png                fig11.png
)
#    ZATOM_bifur_analysis_MLT_S.png                fig10.png
#    ZATOM_bifur_analysis_MLT_S_marks.png          fig11.png

N=$(( ${#name_pairs[@]} / 2 ))
echo "We have $N figure(s) to rename."
for i in $( seq 1 $N ) ; do
    src_file="${name_pairs[$(( (i-1) * 2 + 0 ))]}"
    dst_file="${name_pairs[$(( (i-1) * 2 + 1 ))]}"
    echo "$src_file => $dst_file"
    cp figures/$src_file final_figures/$dst_file 
done

echo "Done."
