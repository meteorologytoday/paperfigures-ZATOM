#!/bin/bash


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data/batch_basin_width/output_redo_tanh/CM_basin_width00000050_pos/lb8  \
   ./data/continuation_data/batch_basin_width/output_redo_tanh/CM_basin_width00000030_pos/lb8  \
    --legend '$\Delta\lambda_\mathrm{tot}=50^\circ$' '$\Delta\lambda_\mathrm{tot}=30^\circ$' \
    --legend-coor 0.1266 8.828 \
                  0.0884 8.119 \
    --F-rng -0.01 0.3 --psi-rng 5 10 --metric max --output-bifur "figures/bifurcation_basin_width.png"
