#!/bin/bash

echo "Plotting regimes..."

julia src/plot_regime_diagrams_comparison.jl  0.0 0.00 "plot_shiftedX" "(a)" &
julia src/plot_regime_diagrams_comparison.jl  0.5 0.00 "plot_shifted"  "(b)" &
julia src/plot_regime_diagrams_comparison.jl  0.0 0.05 "plot_shifted" "(c)" &

wait
