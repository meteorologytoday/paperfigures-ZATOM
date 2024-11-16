#!/bin/bash


res=8


python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000420_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000460_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000480_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000500_pos/lb$res  \
    --legend '$\xi=-4$' '$\xi=-4.2$' '$\xi=-4.4$' '$\xi=-4.6$' '$\xi=-4.8$' '$\xi=-5$' \
    --no-legend \
    --auto-color \
    --param gamma \
    --param-rng -0.01 0.16  \
    --varnames "mode_psi" \
    --mode-psi-rng 2.0 5.0 \
    --output-bifur "figures/ZATOM_dense_gamma-psi.svg" \
    --legend-loc "center right" \
    --residue-threshold 3e-10 \
    --text '$\xi=-4$' '$-4.2$' '$-4.4$' '$-4.6$' '$-4.8$' '$-5$' \
    --text-ax  0 0 0 0 0 0 \
    --text-pos \
               0.125 3.9 \
               0.125 3.40 \
               0.090 2.25 \
               0.065 2.25 \
               0.045 2.25 \
               0.030 2.25 \
    --title '' \
    --put-var-on-yaxis \
    --no-display


