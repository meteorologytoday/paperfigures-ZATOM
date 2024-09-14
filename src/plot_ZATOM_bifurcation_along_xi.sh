#!/bin/bash

res=8

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000005_pos/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000025_pos/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000050_neg/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000075_neg/lb$res  \
   ./data/continuation_data_fixed_gamma_20240908/batch_60/output_redo_balanced_tanh/CM_gamma00000100_neg/lb$res  \
    --legend '$\gamma'"'"'=0.005$' '$\gamma'"'"'=0.025$' '$\gamma'"'"'=0.05$' '$\gamma'"'"'=0.075$' '$\gamma'"'"'=0.1$' \
    --param     xi  \
    --param-rng -7.10 0.1  \
    --psi-rng -2 10 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --residue-threshold 3e-10  \
    --offset-marks 0.01 \
    --output-bifur "figures/ZATOM_bifur_analysis_fixed_gamma.png"


#    --marks 0.0         2.09   \
#            0.1146      3.685  \
#            0.1146     -0.751  \
#            0.1146      0.801  \
#            0.1146     -0.586  \
#    --output-marks "figures/ZATOM_bifur_analysis_xi_marks.png" \

