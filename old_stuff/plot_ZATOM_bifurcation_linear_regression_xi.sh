#!/bin/bash

res=8

python3 src/linear_regression_analysis.py --folder \
   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
    --legend '$\xi'"'"'=1.0$' \
    --colors "black" "orange" "red" "darkgreen" "blue" "purple"  \
    --gamma-rng 0.07 0.13  \
    --output-bifur "figures/ZATOM_bifur_analysis_linear_regression_xi.png" \
    --residue-threshold 3e-10 
    
    #--gamma-rng -0.01 0.16  \

#   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi-0000050_pos/lb$res  \
#   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000000_pos/lb$res  \
#   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000050_pos/lb$res  \
#   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
#   ./data/continuation_data_0531/batch_std20/output_redo_balanced_tanh/CM_xi00000150_pos/lb$res  \
#    --legend '$\xi'"'"'=0.5$' '$\xi'"'"'=0$' '$\xi'"'"'=-0.5$' '$\xi'"'"'=-1$' '$\xi'"'"'=-1.5$' \

