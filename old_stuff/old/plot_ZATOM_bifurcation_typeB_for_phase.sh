#!/bin/bash

res=5

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000075_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000080_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000085_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000090_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000095_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000100_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000105_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000110_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000115_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_balanced_tanh/CM_xi00000200_pos/lb$res  \
    --legend '$\xi=.75$' '$\xi=.80$' '$\xi=.85$' '$\xi=.90$' '$\xi=.95$' '$\xi=1$' '$\xi=1.05$' '$\xi=1.10$' '$\xi=1.15$' '$\xi=2$'\
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.6 --y-rng 0 12 --metric fixed --metric-dep 1000.0 --metric-lat 50 \
    --residue-threshold 3e-10
            
