#!/bin/bash

julia src/plot_etb_dydt.jl --title '(a) $\epsilon = 0$'   --eps 0.0 --output "figures/figure-etb_dydt-a.svg" &
julia src/plot_etb_dydt.jl --title '(b) $\epsilon = 0.5$' --eps 0.5 --output "figures/figure-etb_dydt-b.svg" &


wait
