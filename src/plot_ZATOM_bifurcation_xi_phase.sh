#!/bin/bash

res=8

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000050_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000075_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000080_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000085_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000090_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000095_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000105_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000110_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000115_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000120_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000125_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000130_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000135_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000140_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000145_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000150_pos/lb$res  \
    --legend '0.5' '0.75' '0.80' '0.85' '0.90' '0.95' '1' '1.05' '1.1' '1.15' '1.2' '1.25' '1.3' '1.35' '1.4' '1.45' '1.5' \
    --auto-color \
    --residue-threshold 3e-10  \
    --metric fixed \
    --metric-dep 1000 \
    --metric-lat 60
#    --gamma-rng -0.01 0.2  \
#    --psi-fixed-rng -2 10 \
#    --s1000-rng  2 3   \
#    --db_ew-rng -3.5 4    \
#    --db_ns-rng  3 6    \
#    --cvt_e-rng -0.01 0.2    \
#    --cvt_w-rng -0.01 0.2    \
#    --offset-marks 0.01 \
#    --output-bifur "figures/ZATOM_bifur_analysis_xi_phase.png" \


