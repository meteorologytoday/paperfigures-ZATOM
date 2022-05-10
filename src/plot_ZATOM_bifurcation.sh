#!/bin/bash


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data/batch_std50/output_redo_tanh/CM_xi-0000150_pos/lb8  \
   ./data/continuation_data/batch_std50/output_redo_tanh/CM_xi-0000100_pos/lb8  \
   ./data/continuation_data/batch_std50/output_redo_tanh/CM_xi00000000_pos/lb8  \
   ./data/continuation_data/batch_std50/output_redo_tanh/CM_xi00000100_pos/lb8  \
    --legend '$\xi=-1.5$' '$\xi=-1$' '$\xi=0$' '$\xi=1$' \
    --legend-coor 0.1511 8.990 \
                  0.1504 8.216 \
                  0.1091 7.760 \
                  0.0232 7.809 \
    --F-rng -0.01 0.3 --psi-rng 5 10 --metric max --output-bifur "figures/std50_bifurcation.png" \
    --marks 0.0      9.468 \
            0.2376   7.916 \
            0.2376   6.460 \
    --output-marks "figures/std50_standard_marks.png"
            
