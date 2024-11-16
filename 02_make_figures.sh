#!/bin/bash

source 00_setup.sh
source 98_trapkill.sh

plot_codes=(

    # Figure 3
    $jl $src_dir/30_plot_forcing.jl

    # Figures 4 and 6
    $bs $src_dir/40_plot_ZATOM_bifurcation_diag.sh

    # Figure 5
    $bs $src_dir/50_plot_ZATOM_diff.sh

    # Figure 7
    $bs $src_dir/70_plot_ZATOM_bifurcation_dense_gamma-psi.sh
    $bs $src_dir/71_plot_ZATOM_bifurcation_dense_xi-psi.sh


    # Figure 8 
    $bs $src_dir/80_plot_etb_dydt.sh

    # Figure 9
    $jl $src_dir/90_plot_etb_bifur_p.jl
    $jl $src_dir/91_plot_etb_bifur_xi.jl
    
    # Figure 10
    $jl $src_dir/100_plot_regime_diagrams_comparison.jl

    # Figure D1 
    $jl $src_dir/D10_plot_etb_bifur_phase.jl
)

# Some code to download data and extract them

mkdir figures

nproc=5
N=$(( ${#plot_codes[@]} / 2 ))
echo "We have $N file(s) to run..."
for i in $( seq 1 $(( ${#plot_codes[@]} / 2 )) ) ; do

    {
        PROG="${plot_codes[$(( (i-1) * 2 + 0 ))]}"
        FILE="${plot_codes[$(( (i-1) * 2 + 1 ))]}"

        echo "=====[ Running file: $FILE ]====="
        eval "$PROG $FILE" 
    } &

    proc_cnt=$(( $proc_cnt + 1))

    if (( $proc_cnt >= $nproc )) ; then
        echo "Max proc reached: $nproc"
        wait
        proc_cnt=0
    fi

done


wait

echo "Figures generation is complete."
echo "Please run 03_postprocess_figures.sh to postprocess the figures."
