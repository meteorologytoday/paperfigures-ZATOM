#!/bin/bash

echo "Plotting regimes..."

#julia src/plot_regime_diagrams_comparison.jl  0.0 0.00 "plot_shiftedX" "(a)" &
julia src/plot_regime_diagrams_comparison.jl  0.65 0.00 "plot_shifted"  "(b)" &
#julia src/plot_regime_diagrams_comparison.jl   0.5 0.03 "plot_shifted" "(c)" &

wait
