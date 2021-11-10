#!/bin/bash

### run scripts for all models in parallel

Rscript scripts/PureASM_N0grp_L0mouse.R &
Rscript scripts/PureASM_N0mouse_L0grp.R &
Rscript scripts/PureASM_N0mouse_L0mouse.R &
Rscript scripts/PureASM_N0mouse.R &
wait
