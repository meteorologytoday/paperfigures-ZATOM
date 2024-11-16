#!/bin/bash

res=8

#    --plot-psi-rng 2.0 5.0 \
#    --plot-chi-rng 0.0 2.6 \

python3 src/plot_bifur_analysis_diag.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --plot-param-rng -0.005 0.125 \
    --param gamma \
    --param-rng -0.01 0.12  \
    --param-scale 0.03 \
    --psi-scale 0.2 \
    --mode 1 \
    --label-offset 0.005 \
    --label-idx 17 43 \
    --label-sides R L \
    --labels '$S_{\mathrm{on}}$' '$S_{\mathrm{off}}$'     \
    --output "figures/ZATOM_bifur_diag_xi-440.svg" \
    --nrow 2 \
    --suptitle '$ \xi = -4.4 $' \
    --no-display
 
python3 src/plot_bifur_analysis_diag.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
    --plot-param-rng -0.005 0.125 \
    --param gamma \
    --param-rng -0.01 0.12  \
    --param-scale 0.03 \
    --psi-scale 0.2 \
    --mode 1 \
    --output "figures/ZATOM_bifur_diag_xi-400.svg" \
    --nrow 2 \
    --suptitle '$ \xi = -4 $' \
    --no-display
   
#./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
#    --dont-plot-chi --dont-plot-psi \

