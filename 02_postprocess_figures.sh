#!/bin/bash


echo "Making output directory 'final_figures'..."
mkdir final_figures


echo "Making final figures... "


# Merging two sub-figures
convert \( figures/figure-stommel_bifurcation_phase.png \) \
    \( figures/figure-ZATOM_bifurcation_phase.png \) -gravity center -append \
     figures/merged-phase-diagram.png



name_pairs=(
    bifurcation_xi_typeA.png                    fig03.png
    bifurcation_xi_marks_typeA.png              fig04.png
    bifurcation_xi_typeB.png                    fig05.png
    bifurcation_xi_marks_typeB.png              fig06.png
    figure-reduced_stommel.png                  fig07.png
    figure-stommel_bifurcation_analytical.png   fig08.png
    merged-phase-diagram.png                    fig09.png
)

N=$(( ${#name_pairs[@]} / 2 ))
echo "We have $N figure(s) to rename."
for i in $( seq 1 $(( ${#plot_codes[@]} / 2 )) ) ; do
    src_file="${name_pairs[$(( (i-1) * 2 + 0 ))]}"
    dst_file="${name_pairs[$(( (i-1) * 2 + 1 ))]}"
    echo "$src_file => $dst_file"
    cp figures/$src_file final_figures/$dst_file 
done

echo "Done."
