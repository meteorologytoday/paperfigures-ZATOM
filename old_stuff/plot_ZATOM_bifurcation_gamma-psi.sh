#!/bin/bash

res=8

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --legend '$\xi=-4$' '$\xi=-4.4$' \
    --colors "black" "red" \
    --param gamma \
    --param-rng -0.01 0.16  \
    --mode1-psi-rng 1.3 4.6 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --offset-marks 0.005 \
    --varnames "mode1_psi" "mode1_ZOC" "mode1_dq" "mode1_chi_dbdz" "mode1_ui_adv" \
    --output-bifur "figures/ZATOM_bifur_analysis_xi.svg" \
    --marks  0.0826     2.37   \
             0.04583    1.547  \
    --mark-labels '$S_{R}$' '$S_{L}$'     \
    --mark-sides  "R" "L" \
    --output-marks "figures/ZATOM_bifur_analysis_xi_marks.svg" \
    --residue-threshold 3e-10 \
    --legend-loc "center right" \
    --ncol 5 \
    --const-A 1.35e-11 \
    --const-C 9.57e-19 \
    --const-beta 0.8 \
    --const-tau  $(( 86400 * 864 )) 

