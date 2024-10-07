#!/bin/bash

source 00_setup.sh


finalfig_pdf_dir=$finalfig_dir/pdf
finalfig_png_dir=$finalfig_dir/png
finalfig_svg_dir=$finalfig_dir/svg


echo "Making output directory '${finalfig_dir}'..."
mkdir -p $finalfig_dir
mkdir -p $finalfig_pdf_dir
mkdir -p $finalfig_png_dir
mkdir -p $finalfig_svg_dir


echo "Making final figures... "


cp $fig_static_dir/model.svg $fig_dir/
cp $fig_static_dir/model_physics.svg $fig_dir/




echo "Doing merging : analytical extended two-box model"
svg_stack.py                        \
    --direction=h                   \
    figures/figure-etb_bifur_p.svg  \
    figures/figure-etb_bifur_xi.svg \
    > figures/merged-stommel_bifurcation_analytical.svg


#echo "Doing merging : regime diagrams"
#convert \( figures/regime_diagrams_comparison_xi_0.00_gamma_0.00.png \) \
#    \( figures/regime_diagrams_comparison_xi_0.65_gamma_0.00.png \) -gravity North +append \
#     figures/regime_diagrams_comparison.png


echo "Doing merging : freshwater forcing and cartoon"
svg_stack.py \
    --direction=h \
    $fig_dir/quantitative-forcing.svg \
    $fig_static_dir/fwf_cartoon.svg \
    > $fig_dir/merged-forcing.svg

echo "Doing merging : extended two-box model dydt plots"
svg_stack.py \
    --direction=h \
    figures/figure-etb_dydt-a.svg \
    figures/figure-etb_dydt-b.svg \
    > figures/merged-etb_dydt.svg

echo "Doing merging : freshwater forcing and cartoon"
svg_stack.py \
    --direction=h \
    $fig_dir/ZATOM_dense_gamma-psi.svg \
    $fig_dir/ZATOM_dense_xi-psi.svg \
    > $fig_dir/merged-ZATOM_dense_bifur.svg

name_pairs=(
    model_physics.svg                             fig01
    model.svg                                     fig02
    merged-forcing.svg                            fig03
    ZATOM_bifur_analysis_xi.svg                   fig04
    ZATOM_state_diff.svg                          fig05
    merged-ZATOM_dense_bifur.svg                  fig06
    merged-etb_dydt.svg                           fig07
    merged-stommel_bifurcation_analytical.svg     fig08
    regime_diagrams_comparison.svg                fig09
#    ZATOM_bifur_gamma_xi.png                      fig08

#    figure-etb_bifur_phase.png                    fig10
#    fwf_dijkstra.png                              fig11
#    figure-approx_dijkstra.png                    fig12
)

N=$(( ${#name_pairs[@]} / 2 ))
echo "We have $N figure(s) to rename and convert into pdf files."
for i in $( seq 1 $N ) ; do

    {

    src_file="${name_pairs[$(( (i-1) * 2 + 0 ))]}"
    dst_file_pdf="${name_pairs[$(( (i-1) * 2 + 1 ))]}.pdf"
    dst_file_png="${name_pairs[$(( (i-1) * 2 + 1 ))]}.png"
    dst_file_svg="${name_pairs[$(( (i-1) * 2 + 1 ))]}.svg"
 
    echo "$src_file => $dst_file_svg"
    cp $fig_dir/$src_file $finalfig_svg_dir/$dst_file_svg
   
    echo "$src_file => $dst_file_pdf"
    cairosvg $fig_dir/$src_file -o $finalfig_pdf_dir/$dst_file_pdf

    echo "$src_file => $dst_file_png"
    magick $finalfig_pdf_dir/$dst_file_pdf $finalfig_png_dir/$dst_file_png

    } &
done

wait

echo "Done."


