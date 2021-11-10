#!/bin/bash

### run scripts for all models in parallel

Rscript scripts/APHA_N0mouse_EXP.R &
Rscript scripts/APHA_N0mouse_Sigmoid1.R &
Rscript scripts/APHA_N0mouse_Sigmoid2.R &
wait
