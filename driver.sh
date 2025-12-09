#!/bin/sh --login

# loading modules (this is for GAEA C6)
module use -a /gpfs/f6/bil-fire10-oar/world-shared/gge/Miniforge3/modulefiles
module load Miniforge3/24.11.3-2
module load pyDAmonitor/1.0.0


module list


# horizontal_mpasjedi_diff
mkdir -p fig_horizontal_mpasjedi_diff
python ./scripts/horizontal_mpasjedi_diff.py
