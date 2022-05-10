#!/bin/bash

py=python3
jl=julia
bs=bash
srcdir=src

plot_codes=(
    $jl $srcdir/plot_reduced_stommel.jl
    $jl $srcdir/plot_bifurcation_analytical.jl
    $bs $srcdir/plot_ZATOM_bifurcation.sh
)

# Some code to download data and extract them


mkdir figures

N=$(( ${#plot_codes[@]}  ))
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
