#!/bin/bash


python3 src/plot_scans_and_snapshots.py --folder \
   ./data/continuation_data/batch_cva50/output_redo_tanh/CM_cva_delta00000005_pos/lb8  \
   ./data/continuation_data/batch_cva50/output_redo_tanh/CM_cva_delta00000010_pos/lb8  \
   ./data/continuation_data/batch_cva50/output_redo_tanh/CM_cva_delta00000020_pos/lb8  \
    --legend '$\Delta_\mathrm{cm}=5\times10^{-5}\mathrm{m}/\mathrm{s}^2$' '$\Delta_\mathrm{cm}=1\times10^{-4}\mathrm{m}/\mathrm{s}^2$' '$\Delta_\mathrm{cm}=2\times10^{-4}\mathrm{m}/\mathrm{s}^2$' \
    --legend-coor 0.1761 9.410 \
                  0.2006 8.575 \
                  0.0489 7.702 \
    --F-rng -0.01 0.3 --psi-rng 5 10 --metric max --output-bifur "figures/bifurcation_cva.png" 
            
