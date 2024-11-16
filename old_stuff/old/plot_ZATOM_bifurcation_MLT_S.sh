#!/bin/bash

res=8

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_neg/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
    --legend '$H_S=58\mathrm{m}$' '$H_S=204\mathrm{m}$' '' '$H_S=400\mathrm{m}$' \
    --colors "darkorange" "darkgreen" "darkgreen" "blue" \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 1.5 --y-rng -2 10 --metric fixed --metric-dep 1000.0 --metric-lat 60 \
    --output-bifur "figures/bifurcation_xi_typeA.png" \
    --marks 0.0        5.360 \
            0.3590     1.360 \
            0.3590     -1.61 \
            0.8508      5.47 \
            0.8508     -0.234 \
    --marks-pos "right" "right" "bottom" "top" "top" \
    --output-marks "figures/ZATOM_bifurcation_MLT_S.png" \
    --residue-threshold 3e-10 &


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_neg/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
    --legend '$H_S=58\mathrm{m}$' '$H_S=204\mathrm{m}$' '' '$H_S=400\mathrm{m}$' \
    --colors "darkorange" "darkgreen" "darkgreen" "blue" \
    --x-var "gamma" --y-var "db_ns" \
    --x-rng -0.01 1.5 --y-rng NaN NaN --metric max \
    --residue-threshold 3e-10 &


 python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_neg/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
    --legend '$H_S=58\mathrm{m}$' '$H_S=204\mathrm{m}$' '' '$H_S=400\mathrm{m}$' \
    --colors "darkorange" "darkgreen" "darkgreen" "blue" \
    --x-var "gamma" --y-var "db_ew" \
    --x-rng -0.01 1.5 --y-rng NaN NaN --metric max \
    --residue-threshold 3e-10 &


 
python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_neg/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
    --legend '$H_S=58\mathrm{m}$' '$H_S=204\mathrm{m}$' '' '$H_S=400\mathrm{m}$' \
    --colors "darkorange" "darkgreen" "darkgreen" "blue" \
    --x-var "gamma" --y-var "s1000" \
    --x-rng -0.01 1.5 --y-rng NaN NaN --metric max \
    --residue-threshold 3e-10 &

wait
 
            
