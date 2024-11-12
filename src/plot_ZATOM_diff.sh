#!/bin/bash

res=8

python3 src/plot_diff_states.py --folder \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
   ./data/continuation_data_fixed_xi_20240908/batch_60/output_redo_balanced_tanh/CM_xi-0000440_pos/lb$res  \
    --diff-idx 17 43 \
    --titles '$S_{\mathrm{on}}$' '$\Delta S = S_{\mathrm{off}} - S_{\mathrm{on}}$' \
    --output "figures/ZATOM_state_diff.svg" \
    --no-display

