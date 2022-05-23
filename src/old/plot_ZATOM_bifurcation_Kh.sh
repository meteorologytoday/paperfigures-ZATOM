#!/bin/bash


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data/batch_Kh/output_redo_tanh/CM_Kh02000000_pos/lb8  \
   ./data/continuation_data/batch_Kh/output_redo_tanh/CM_Kh03000000_pos/lb8  \
   ./data/continuation_data/batch_Kh/output_redo_tanh/CM_Kh04000000_pos/lb8  \
   ./data/continuation_data/batch_Kh/output_redo_tanh/CM_Kh05000000_pos/lb8  \
    --legend '$K_H=2\times10^4 \mathrm{m}^2/\mathrm{s}$' '$K_H=3\times10^4 \mathrm{m}^2/\mathrm{s}$' '$K_H=4\times10^4 \mathrm{m}^2/\mathrm{s}$' '$K_H=5\times10^4 \mathrm{m}^2/\mathrm{s}$' \
    --legend-coor 0.1628 11.77 \
                  0.0870 10.38 \
                  0.1543  8.85 \
                  0.0800  7.43 \
    --colors "blue" "darkorange" "red" "darkgreen" \
    --F-rng -0.01 0.36 --psi-rng 5 12.5 --metric max --output-bifur "figures/bifurcation_Kh.png"
            
