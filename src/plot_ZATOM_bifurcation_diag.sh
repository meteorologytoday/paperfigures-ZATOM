#!/bin/bash

res=8

python3 src/plot_bifur_analysis_diag.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --param gamma \
    --param-rng -0.01 0.16  \
    --param-scale 0.03 \
    --psi-scale 0.2 \
    --mode 1 \
    --label-offset 0.005 \
    --label-idx 17 43 \
    --label-sides R L \
    --labels '$S_{A}$' '$S_{B}$'     \
    --output "figures/ZATOM_bifur_diag_xi-440.svg" \
    --no-display
 
python3 src/plot_bifur_analysis_diag.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
    --param gamma \
    --param-rng -0.01 0.16  \
    --param-scale 0.03 \
    --psi-scale 0.2 \
    --mode 1 \
    --output "figures/ZATOM_bifur_diag_xi-400.svg" \
    --no-display
   
#./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
#    --dont-plot-chi --dont-plot-psi \

