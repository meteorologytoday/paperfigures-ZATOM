#!/bin/bash


echo "Making output directory 'final_figures'..."
mkdir final_figures


echo "Making final figures... "


# Merging two sub-figures
#convert \( figures/figure-stommel_bifurcation_phase.png \) \
#    \( figures/figure-ZATOM_bifurcation_phase.png \) -gravity center -append \
#     figures/merged-phase-diagram.png

name_pairs=(
    figure-forcing.png                            fig03.png
    figure-reduced_stommel.png                    fig04.png
    figure-stommel_bifurcation_analytical.png     fig05.png
    figure-stommel_bifurcation_phase.png          fig06.png
    ZATOM_bifur_analysis_xi.png                   fig07.png
    ZATOM_bifur_analysis_xi_marks.png             fig08.png
    figure-ZATOM_bifurcation_phase.png            fig09.png
    ZATOM_bifur_analysis_MLT_S.png                fig10.png
    ZATOM_bifur_analysis_MLT_S_marks.png          fig11.png
)

N=$(( ${#name_pairs[@]} / 2 ))
echo "We have $N figure(s) to rename."
for i in $( seq 1 $N ) ; do
    src_file="${name_pairs[$(( (i-1) * 2 + 0 ))]}"
    dst_file="${name_pairs[$(( (i-1) * 2 + 1 ))]}"
    echo "$src_file => $dst_file"
    cp figures/$src_file final_figures/$dst_file 
done

echo "Done."
