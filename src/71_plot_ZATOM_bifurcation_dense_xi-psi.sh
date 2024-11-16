#!/bin/bash

res=8

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000005_pos/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000025_neg/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000050_neg/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000075_neg/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000100_neg/lb$res  \
    --no-legend                \
    --text \
        '$\gamma=.005\,\mathrm{Sv}$'   \
        '$.025\,\mathrm{Sv}$'   \
        '$.050\,\mathrm{Sv}$'   \
        '$.075\,\mathrm{Sv}$'   \
        '$.100\,\mathrm{Sv}$'   \
    --text-ax 0 0 0 0 0  \
    --text-pos   \
        -1 2.90    \
        -1 3.45    \
        -1 3.85    \
        -2 4.0    \
        -3.4 4.15    \
    --auto-color            \
    --param xi              \
    --param-rng -7.5  0     \
    --varnames "mode_psi"  \
    --mode-psi-rng 2.0 5.0 \
    --output-bifur "figures/ZATOM_dense_xi-psi.svg" \
    --thumbnail-skip 1 \
    --legend-loc "upper left" \
    --residue-threshold 3e-10 \
    --title "" \
    --put-var-on-yaxis \
    --no-display

