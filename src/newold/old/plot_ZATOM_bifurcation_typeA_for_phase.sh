#!/bin/bash

res=5

python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000100_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000075_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000050_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi-0000025_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000000_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000025_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000050_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000075_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000100_pos/lb$res  \
   ./data/continuation_data_0522/batch_std20/output_redo_raw_tanh/CM_xi00000200_pos/lb$res  \
    --legend '$\xi=-1$' '$\xi=-.75$' '$\xi=-.5$' '$\xi=-.25$' '$\xi=0$' '$\xi=.25$' '$\xi=.5$' '$\xi=.75$' '$\xi=1$' '$\xi=2$' \
    --x-var "gamma" --y-var "psi" \
    --x-rng -0.01 0.6 --y-rng -8 6 --metric fixed --metric-dep 1000.0 --metric-lat 50 \
    --residue-threshold 3e-10
            
