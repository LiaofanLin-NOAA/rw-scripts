#!/bin/bash
#SBATCH -M c6
#SBATCH -A bil-pmp 
#SBATCH -p batch
#SBATCH -t 00:30:00 
#SBATCH -J multivar-inc
#SBATCH -N 1
#SBATCH -n 1





# Set up for the script (user's specification)
# --------------------------------------------
# python path
export PATH_PYDAMONITOR="/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/260801-rrfs-workflow-ens/rrfs-workflow/workflow/sideload/pyDAmonitor" 

# file paths 
# - File A for analysis; File B for background
PATH_RETROA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/260801-rrfs-workflow-ens/OPSROOT/hrly_12km04/stmp/20240527/rrfs_prep_ic_02_v2.1.4/enkf/mem001'
PATH_RETROB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/260801-rrfs-workflow-ens/OPSROOT/hrly_12km04/stmp/20240527/rrfs_fcst_01_v2.1.4/enkf/mem001'

FILEA="$PATH_RETROA/mpasout.nc" 
FILEB="$PATH_RETROB/mpasout.2024-05-27_02.00.00.nc" 

STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/260801-rrfs-workflow-ens/rrfs-workflow/fix/conus12km/conus12km.invariant.nc_L60_GEFS'


# Variables and levels to be processed
# ------------------------------------
# NOTE:
#   Atmospheric variables: level 1 = bottom model level/layer.
#   Land-surface variables: level 1 = top soil level/layer.
VAR_LEVELS=(
    "smois:1"
    "tslb:1"
    "theta:1"
    "qv:1"
)


# Load the python env
# -------------------
source "$PATH_PYDAMONITOR/ush/load_pyDAmonitor.sh"
module list

mkdir -p figures



# Run the script
# --------------
for item in "${VAR_LEVELS[@]}"; do
    var="${item%%:*}"
    level="${item##*:}"

	python mpasjedi-multivar-inc-comparison.py \
	    --filea "$FILEA" \
	    --fileb "$FILEB" \
	    --static "$STATIC" \
		--variable "${var}" \
		--level "${level}"

done


