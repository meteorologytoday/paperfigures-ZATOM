#!/bin/bash

res=8

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000058_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_pos/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000204_neg/lb$res  \
   ./data/continuation_data_0531/batch_MLT_S/output_redo_balanced_tanh/CM_MLT_S00000400_pos/lb$res  \
    --legend '$H_S=58\mathrm{m}$' '$H_S=204\mathrm{m}$' '' '$H_S=400\mathrm{m}$' \
    --colors "red" "darkgreen" "darkgreen" "blue" \
    --gamma-rng -0.01 1.5  \
    --psi-rng -2 12 \
    --s1000-rng 0.5 3.5   \
    --db_ns-rng -2 5.5    \
    --db_ew-rng -5.5 5    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --output-bifur "figures/ZATOM_bifur_analysis_MLT_S.png" \
    --marks 0.3590     1.360 \
            0.3590     -1.61 \
            0.8508      5.47 \
            0.8508     -0.234 \
    --offset-marks 0.05 \
    --output-marks "figures/ZATOM_bifur_analysis_MLT_S_marks.png" \
    --residue-threshold 3e-10 

