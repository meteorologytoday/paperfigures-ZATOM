#!/bin/bash


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi-0000100_pos/lb5  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi-0000050_pos/lb5  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000000_pos/lb5  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000050_pos/lb5  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb5  \
    --legend '$\xi=-1$' '$\xi=-0.5$' '$\xi=0$' '$\xi=0.5$' '$\xi=1$' \
    --colors "darkorange" "red" "darkgreen" "blue" "purple" \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.6 --y-rng 0 12 --metric fixed --metric-dep 1000.0 --metric-lat 50 \
    --output-bifur "figures/bifurcation_xi_typeB.png" \
    --marks 0.433    9.972 \
    --output-marks "figures/bifurcation_xi_marks_typeB.png" \
    --residue-threshold 5e-10
            
