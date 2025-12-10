#!/bin/bash --login

# loading modules (this is for GAEA C6)
module use -a /gpfs/f6/bil-fire10-oar/world-shared/gge/Miniforge3/modulefiles
module load Miniforge3/24.11.3-2
module load pyDAmonitor/1.0.0


module list


# horizontal_mpasjedi_diff
FILEA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km03/com/rrfs/v2.1.2/rrfs.20240502/11/fcst/det/mpasout.2024-05-02_12.00.00.nc' #1h forecast, init at 11utc,    sfc continuous cycling
FILEB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km02/com/rrfs/v2.1.2/rrfs.20240502/11/fcst/det/mpasout.2024-05-02_12.00.00.nc' #1h forecast, init at 11utc, NO sfc continuous cycling
STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L65_GFS'

mkdir -p fig_horizontal_mpasjedi_diff
python ./scripts/horizontal_mpasjedi_diff.py \
    --filea  "$FILEA" \
    --fileb  "$FILEB" \
    --static "$STATIC"


# horizontal_mpasjedi_diff
FILEA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km03/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc' #1h forecast, init at 11utc,    sfc continuous cycling
FILEB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km02/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc' #1h forecast, init at 11utc, NO sfc continuous cycling
STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L65_GFS'

mkdir -p fig_horizontal_mpasjedi_diff
python ./scripts/horizontal_mpasjedi_diff.py \
    --filea  "$FILEA" \
    --fileb  "$FILEB" \
    --static "$STATIC"

