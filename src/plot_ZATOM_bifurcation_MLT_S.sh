#!/bin/bash

res=8

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000125_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000295_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000522_pos/lb$res  \
    --legend '$\xi=-1$' '$\xi=-0.5$' '$\xi=0$' '$\xi=0.5$' '$\xi=1$' \
    --auto-color \
    --colors "darkorange" "red" "darkgreen" "blue" "purple" \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.6 --y-rng -8 6 --metric fixed --metric-dep 1000.0 --metric-lat 50 \
    --output-bifur "figures/bifurcation_xi_typeA.png" \
    --marks 0.0        5.360 \
            0.4551    -2.604 \
            0.4551    -5.436 \
            0.4551     4.943 \
    --marks-pos "right" "right" "bottom" "top"\
    --output-marks "figures/ZATOM_bifurcation_MLT_S.png" \
    --residue-threshold 3e-10
            
