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

cp manual_figures/* figures/

echo "Doing merging : analytical extended two-box model"
convert \( figures/figure-etb_bifur_p.png \) \
    \( figures/figure-etb_bifur_xi.png \) -gravity North +append \
     figures/merged-stommel_bifurcation_analytical.png


#echo "Doing merging : regime diagrams"
#convert \( figures/regime_diagrams_comparison_xi_0.00_gamma_0.00.png \) \
#    \( figures/regime_diagrams_comparison_xi_0.65_gamma_0.00.png \) -gravity North +append \
#     figures/regime_diagrams_comparison.png


echo "Doing merging : freshwater forcing and cartoon"
convert \( figures/figure-forcing.png -gravity North -background white -splice 0x200 -pointsize 120 -annotate +0+120 '(a)' \) \
    \( figures/fwf_zatom.png -resize 70% -background white -splice 0x200 -pointsize 120 -annotate +0+120 '(b)' \) \
     -gravity North +append \
     figures/merged-forcing.png

echo "Doing merging : extended two-box model dydt plots"
convert figures/figure-etb_dydt-a.png figures/figure-etb_dydt-b.png -gravity North +append \
     figures/merged-etb_dydt.png


name_pairs=(
    MOC_dynamics.png                              fig01.png
    ZATOM_design.png                              fig02.png 
    merged-forcing.png                            fig03.png
    ZATOM_bifur_analysis_xi.png                   fig04.png
    ZATOM_bifur_analysis_xi_marks.png             fig05.png
    merged-etb_dydt.png                           fig06.png
    merged-stommel_bifurcation_analytical.png     fig07.png
    ZATOM_bifur_gamma_xi.png                      fig08.png
    regime_diagrams_comparison.png                fig09.png
    figure-etb_bifur_phase.png                    fig10.png
    fwf_dijkstra.png                              fig11.png
    figure-approx_dijkstra.png                    fig12.png
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
