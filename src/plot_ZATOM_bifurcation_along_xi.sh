#!/bin/bash

res=8

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000000_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000002_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000004_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000006_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000008_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000010_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000012_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000014_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000016_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000018_neg/lb$res  \
   ./data/continuation_data_xi_072422/batch_std/output_redo_balanced_tanh/CM_gamma00000020_neg/lb$res  \
    --legend '$\gamma'"'"'=0$' '$\gamma'"'"'=0.02$' '$\gamma'"'"'=0.04$' '$\gamma'"'"'=0.06$' '$\gamma'"'"'=0.08$' '$\gamma'"'"'=0.1$' '$\gamma'"'"'=0.12$' '$\gamma'"'"'=0.14$' '$\gamma'"'"'=0.16$' '$\gamma'"'"'=0.18$' '$\gamma'"'"'=0.2$' \
    --param     xi  \
    --param-rng -2.00 1.01  \
    --psi-rng -2 10 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --residue-threshold 3e-10  \
    --offset-marks 0.01 \
    --output-bifur "figures/ZATOM_bifur_analysis_xi.png"


#    --marks 0.0         2.09   \
#            0.1146      3.685  \
#            0.1146     -0.751  \
#            0.1146      0.801  \
#            0.1146     -0.586  \
#    --output-marks "figures/ZATOM_bifur_analysis_xi_marks.png" \

