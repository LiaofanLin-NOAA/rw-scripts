#!/bin/bash
#SBATCH -M c6
#SBATCH -A bil-pmp 
#SBATCH -p batch
#SBATCH -t 00:30:00 
#SBATCH -J multivar-inc
#SBATCH -N 1
#SBATCH -n 1


# Load the python env
# -------------------
export PATH_PYDAMONITOR="/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/260801-rrfs-workflow-ens/rrfs-workflow/workflow/sideload/pyDAmonitor"

source $PATH_PYDAMONITOR/ush/load_pyDAmonitor.sh
module list


# Set up for the script
# ---------------------
PATH_RETROA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km03'
PATH_RETROB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km02'

FILEA="$PATH_RETROA/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc" 
FILEB="$PATH_RETROB/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc" 

STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L65_GFS'



# Run the script
# --------------
python mpasjedi-multivar-inc-comparison.py \
    --filea "$FILEA" \
    --fileb "$FILEB" \
    --static "$STATIC"