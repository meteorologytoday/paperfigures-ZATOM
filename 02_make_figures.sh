#!/bin/bash

source 00_setup.sh
source 98_trapkill.sh

plot_codes=(
    $jl $src_dir/plot_forcing.jl
#    $bs $src_dir/plot_ZATOM_bifurcation_gamma-psi.sh
    $bs $src_dir/plot_ZATOM_bifurcation_diag.sh
    $bs $src_dir/plot_ZATOM_bifurcation_dense_xi-psi.sh
    $bs $src_dir/plot_ZATOM_bifurcation_dense_gamma-psi.sh
    $bs $src_dir/plot_ZATOM_diff.sh
    $bs $src_dir/plot_etb_dydt.sh
    $jl $src_dir/plot_etb_bifur_p.jl
    $jl $src_dir/plot_etb_bifur_xi.jl
    $jl $src_dir/plot_etb_bifur_phase.jl
    $jl $src_dir/plot_regime_diagrams_comparison.jl
#    $bs $src_dir/plot_ZATOM_bifurcation_MLT_S.sh
#    $py $src_dir/plot_ZATOM_bifurcation_gamma_xi.py
    $jl $src_dir/plot_approx_dijkstra.jl
)

#plot_codes=(
#    $bs $src_dir/plot_ZATOM_bifurcation_xi.sh
#)


#$jl $src_dir/plot_ZATOM_bifurcation_phase.jl
#$jl $src_dir/plot_ZATOM_bifurcation_phase_HS.jl

# Some code to download data and extract them


mkdir figures

nproc=2
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
