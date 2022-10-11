#!/bin/bash

py=python3
jl=julia
bs=bash
srcdir=src

plot_codes=(
    $jl $srcdir/plot_forcing.jl
    $bs $srcdir/plot_etb_dydt.sh
    $jl $srcdir/plot_etb_bifur_p.jl
    $jl $srcdir/plot_etb_bifur_xi.jl
    $jl $srcdir/plot_etb_bifur_phase.jl
    $jl $srcdir/plot_regime_diagrams_comparison.jl
    $bs $srcdir/plot_ZATOM_bifurcation_xi.sh
    $bs $srcdir/plot_ZATOM_bifurcation_MLT_S.sh
    $py $srcdir/plot_ZATOM_bifurcation_gamma_xi.py
    $jl $srcdir/plot_approx_dijkstra.jl
)

#plot_codes=(
#    $bs $srcdir/plot_ZATOM_bifurcation_xi.sh
#)


#$jl $srcdir/plot_ZATOM_bifurcation_phase.jl
#$jl $srcdir/plot_ZATOM_bifurcation_phase_HS.jl

# Some code to download data and extract them


mkdir figures

N=$(( ${#plot_codes[@]} / 2 ))
echo "We have $N file(s) to run..."
for i in $( seq 1 $(( ${#plot_codes[@]} / 2 )) ) ; do
    PROG="${plot_codes[$(( (i-1) * 2 + 0 ))]}"
    FILE="${plot_codes[$(( (i-1) * 2 + 1 ))]}"
    echo "=====[ Running file: $FILE ]====="
    eval "$PROG $FILE" &
done


wait

echo "Figures generation is complete."
echo "Please run 02_postprocess_figures.sh to postprocess the figures."
