#!/bin/bash

res=8

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi-0000050_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000000_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000050_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000090_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000125_pos/lb$res  \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000150_pos/lb$res  \
    --legend '-0.5' '0' '0.5' '.90' '1' '1.25' '1.5' \
    --auto-color \
    --colors "darkorange" "red" "darkgreen" "blue" "purple" \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.20 --y-rng -1.5 6 --metric fixed --metric-dep 1000.0 --metric-lat 60 \
    --output-bifur "figures/bifurcation_xi_typeA.png" \
    --marks 0.0         2.09   \
            0.1146      3.685  \
            0.1146      0.801  \
            0.1146     -0.586  \
    --marks-pos "right" "right" "bottom" "top"  \
    --output-marks "figures/ZATOM_bifurcation_xi.png" \
    --residue-threshold 3e-10
            
