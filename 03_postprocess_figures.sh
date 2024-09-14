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
convert \( figures/figure-etb_bifur_p.png \) \
    \( figures/figure-etb_bifur_xi.png \) -gravity North +append \
     figures/merged-stommel_bifurcation_analytical.png


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
convert figures/figure-etb_dydt-a.png figures/figure-etb_dydt-b.png -gravity North +append \
     figures/merged-etb_dydt.png


if [ ] ; then
echo "Figure 7: Merge Fourier analysis... "
fixed_dSST=100
fixed_wnm=010
svg_stack.py \
    --direction=v \
    $fig_dir/spectral_analysis_linearity_on_dSST/linearity_on_dSST_lab_FIXEDDOMAIN_SST_sine_WETLWSW_wnm${fixed_wnm}_MYNN25_hr120-240.svg \
    $fig_dir/spectral_analysis_tracking_wnm/spectral_analysis_lab_FIXEDDOMAIN_SST_sine_WETLWSW_dT${fixed_dSST}_MYNN25_hr120-240.svg \
    > $fig_dir/merged-spectral_analysis_wnm${fixed_wnm}_dSST${fixed_dSST}_MYNN25_hr120-240.svg



echo "Figure 2: Merge experiment design and vertical profile..."
svg_stack.py \
    --direction=h \
    $fig_static_dir/experiment_design_3.svg \
    $fig_dir/input_sounding_woML.svg \
    > $fig_dir/merged-exp.svg

echo "Figure 3: Merge snapshots... "

for dT in 100 300; do
for wnm in 004 010 ; do
    
    svg_stack.py \
        --direction=h \
        $fig_dir/snapshots_dhr-120/lab_FIXEDDOMAIN_SST_sine_DRY/case_mph-off_wnm${wnm}_U20_dT${dT}_MYNN25/snapshot-part1_120-240.svg \
        $fig_dir/snapshots_dhr-120/lab_FIXEDDOMAIN_SST_sine_WETLWSW/case_mph-on_wnm${wnm}_U20_dT${dT}_MYNN25/snapshot-part1_120-240.svg \
        > $fig_dir/merged-snapshot_wnm${wnm}_U20_dT${dT}_part1.svg

    svg_stack.py \
        --direction=v \
        $fig_dir/snapshots_dhr-120/lab_FIXEDDOMAIN_SST_sine_DRY/case_mph-off_wnm${wnm}_U20_dT${dT}_MYNN25/snapshot-part2_120-240.svg \
        $fig_dir/snapshots_dhr-120/lab_FIXEDDOMAIN_SST_sine_WETLWSW/case_mph-on_wnm${wnm}_U20_dT${dT}_MYNN25/snapshot-part2_120-240.svg \
        > $fig_dir/merged-snapshot_wnm${wnm}_U20_dT${dT}_part2.svg

done
done



# Merging the phase diagram
echo "Figure 6: Merge phase diagram... "
fixed_dSST=100
fixed_wnm=010
svg_stack.py \
    --direction=v \
    $fig_dir/dF_flux_decomposition_varying_dSST/lab_FIXEDDOMAIN_SST_sine_WETLWSW/dF_flux_decomposition_onefig_wnm${fixed_wnm}_varying_dSST_MYNN25_hr120-240.svg \
    $fig_dir/dF_flux_decomposition_varying_wnm/lab_FIXEDDOMAIN_SST_sine_WETLWSW/dF_flux_decomposition_onefig_dSST${fixed_dSST}_varying_wnm_MYNN25_hr120-240.svg \
    > $fig_dir/merged-dF_flux_decomposition_wnm${fixed_wnm}_dSST${fixed_dSST}_MYNN25_hr120-240.svg

svg_stack.py \
    --direction=v \
    $fig_dir/spectral_analysis_tracking_wnm1/spectral_analysis_lab_FIXEDDOMAIN_SST_sine_WETLWSW_MYNN25_hr120-240.svg \
    $fig_dir/spectral_analysis_tracking_wnm1/spectral_analysis_lab_FIXEDDOMAIN_SST_sine_DRY_MYNN25_hr120-240.svg \
    > $fig_dir/merged-spectral_analysis_tracking_wnm1_MYNN25_hr120-240.svg




echo "Figure 9: Merging the misc phase diagram..."
svg_stack.py \
    --direction=h \
    $fig_dir/phase_misc/lab_FIXEDDOMAIN_SST_sine_WETLWSW/phase_misc_wnm010_varying_dSST_hr120-240.svg \
    $fig_dir/phase_misc/lab_FIXEDDOMAIN_SST_sine_WETLWSW/phase_misc_dSST100_varying_wnm_hr120-240.svg \
    > $fig_dir/merged-phase_misc_hr120-240.svg

fi


name_pairs=(
    model.svg                                     fig01
    model_physics.svg                             fig02 
    merged-forcing.svg                            fig03
#    ZATOM_bifur_analysis_xi.png                   fig04
#    ZATOM_bifur_analysis_xi_marks.png             fig05
#    merged-etb_dydt.png                           fig06
#    merged-stommel_bifurcation_analytical.png     fig07
#    ZATOM_bifur_gamma_xi.png                      fig08
#    regime_diagrams_comparison.png                fig09
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


