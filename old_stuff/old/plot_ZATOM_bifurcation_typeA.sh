#!/bin/bash

res=5

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000100_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000050_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000000_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000050_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000100_pos/lb$res  \
    --legend '$\xi=-1$' '$\xi=-0.5$' '$\xi=0$' '$\xi=0.5$' '$\xi=1$' \
    --colors "darkorange" "red" "darkgreen" "blue" "purple" \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.6 --y-rng -8 6 --metric fixed --metric-dep 1000.0 --metric-lat 50 \
    --output-bifur "figures/bifurcation_xi_typeA.png" \
    --marks 0.0        5.360 \
            0.4551    -2.604 \
            0.4551    -5.436 \
            0.4551     4.943 \
    --marks-pos "right" "right" "bottom" "top"\
    --output-marks "figures/bifurcation_xi_marks_typeA.png" \
    --residue-threshold 3e-10
            
