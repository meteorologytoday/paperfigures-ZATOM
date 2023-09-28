#!/bin/bash

res=8

if [ ] ; then
python3 src/plot_bifur_analysis_new.py --folder \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
    --legend '$\xi=-1$' \
    --colors "black" "orange" "red" "darkgreen" "blue" "purple"  \
    --param gamma \
    --param-rng -0.01 0.2  \
    --psi-fixed-rng -2 10 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.7    \
    --cvt_w-rng -0.01 0.7    \
    --offset-marks 0.01 \
    --output-bifur "figures/ZATOM_bifur_analysis_xi.png" \
    --marks 0.0         2.09   \
            0.1146      3.685  \
            0.1146     -0.751  \
            0.1146      0.801  \
            0.1146     -0.586  \
    --output-marks "figures/ZATOM_bifur_analysis_xi_marks.png" \
    --residue-threshold 3e-10 



python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi-0000100_pos/lb$res  \
    --legend '$\xi=-1$' \
    --colors "darkgreen" \
    --param gamma \
    --param-rng -0.01 0.16  \
    --psi-rng -2 10 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --offset-marks 0.01 \
    --output-bifur "figures/test_ZATOM_bifur_analysis_xi.png" \
    --marks 0.0         2.42   \
            0.1147      4.20  \
            0.1147      1.63  \
            0.1147      1.14  \
    --output-marks "figures/test_ZATOM_bifur_analysis_xi_marks.png" \
    --residue-threshold 3e-10

fi

python3 src/plot_bifur_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi00000050_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi00000000_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi-0000050_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi-0000100_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20230915/batch_std20/output_redo_balanced_tanh/CM_xi-0000150_pos/lb$res  \
    --legend '$\xi=0.5$' '$\xi=0$' '$\xi=-0.5$' '$\xi=-1$' '$\xi=-1.5$' \
    --colors "black" "orange" "red" "darkgreen" "blue" "purple"  \
    --param gamma \
    --param-rng -0.01 0.16  \
    --psi-rng -2 10 \
    --s1000-rng  2 3   \
    --db_ew-rng -3.5 4    \
    --db_ns-rng  3 6    \
    --cvt_e-rng -0.01 0.2    \
    --cvt_w-rng -0.01 0.2    \
    --offset-marks 0.01 \
    --output-bifur "figures/ZATOM_bifur_analysis_xi.png" \
    --marks 0.0         2.42   \
            0.1147      4.20  \
            0.1147      1.63  \
            0.1147      1.14  \
    --output-marks "figures/ZATOM_bifur_analysis_xi_marks.png" \
    --residue-threshold 3e-10 
