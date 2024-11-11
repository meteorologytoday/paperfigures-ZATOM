#!/bin/bash

res=8

if [ ] ; then
python3 src/linear_regression_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --beta-rng 0.0 1.0 \
    --tau-rng 5 $(( 360 * 10 )) \
    --gamma-rng 0.025 0.1 \
    --gamma-scale 0.03 \
    --psi-scale 0.2 \
    --sampling-spacing 0.1 \
    --const-A 1.35e-11 \
    --const-C 9.57e-19 \
    --N-tau  720 \
    --N-beta 11 \
    --output gendata/linear_regression.nc \
    --residue-threshold 3e-10 \
    --no-display
fi

python3 src/linear_regression_analysis.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000400_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --beta-rng 0.0 1.0 \
    --tau-rng $(( 360 * 2 )) $(( 360 * 3 )) \
    --gamma-rng 0.025 0.1 \
    --gamma-scale 0.03 \
    --psi-scale 0.2 \
    --sampling-spacing 0.1 \
    --const-A 1.35e-11 \
    --const-C 9.57e-19 \
    --N-tau  181 \
    --N-beta 21 \
    --output gendata/linear_regression_finer.nc \
    --residue-threshold 3e-10 \
    --no-display

