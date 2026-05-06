#!/bin/bash --login


# =================
# Note and input
# =================
# run "bash driver-post.sh"

# ursa or gaea
hpc_system="ursa"

# =================
# Processes on Ursa
# =================
if [ "${hpc_system}" = "ursa" ]; then

	# Load modules
	module use -a /scratch3/BMC/wrfruc/gge/Miniforge3/modulefiles
	module load Miniforge3/24.11.3-2
	module load pyDAmonitor/1.0.0

	# horizontal_mpasjedi_diff
	FILEA='/scratch4/BMC/wrfruc/llin/ursa/250801-rrfs-workflow/OPSROOT/hrly_12km01/com/rrfs/v2.1.1/rrfs.20240527/02/fcst/det/mpasout.2024-05-27_03.00.00.nc' #1h forecast, init at 11utc,    sfc continuous cycling
	FILEB='/scratch4/BMC/wrfruc/llin/ursa/250801-rrfs-workflow/OPSROOT/hrly_12km01/com/rrfs/v2.1.1/rrfs.20240527/01/fcst/det/mpasout.2024-05-27_02.00.00.nc' #1h forecast, init at 11utc, NO sfc continuous cycling
	STATIC='/scratch4/BMC/wrfruc/llin/ursa/250801-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L60_GFS'

	mkdir -p fig_horizontal_mpasjedi_diff
	python ./scripts-post/horizontal_mpasjedi_diff.py \
	    --filea  "$FILEA" \
	    --fileb  "$FILEB" \
	    --static "$STATIC"
	
	
	
	
fi

# =================
# Processes on Gaea
# =================
if [ "${hpc_system}" = "gaea" ]; then

	# Load modules
	module use -a /gpfs/f6/bil-fire10-oar/world-shared/gge/Miniforge3/modulefiles
	module load Miniforge3/24.11.3-2
	module load pyDAmonitor/1.0.0

	# horizontal_mpasjedi_diff
	FILEA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km03/com/rrfs/v2.1.2/rrfs.20240502/11/fcst/det/mpasout.2024-05-02_12.00.00.nc' #1h forecast, init at 11utc,    sfc continuous cycling
	FILEB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km02/com/rrfs/v2.1.2/rrfs.20240502/11/fcst/det/mpasout.2024-05-02_12.00.00.nc' #1h forecast, init at 11utc, NO sfc continuous cycling
	STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L65_GFS'

	mkdir -p fig_horizontal_mpasjedi_diff
	python ./scripts-post/horizontal_mpasjedi_diff.py \
	    --filea  "$FILEA" \
	    --fileb  "$FILEB" \
	    --static "$STATIC"


	# horizontal_mpasjedi_diff
	FILEA='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km03/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc' #1h forecast, init at 11utc,    sfc continuous cycling
	FILEB='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/OPSROOT/hrly_12km02/com/rrfs/v2.1.2/rrfs.20240502/12/fcst/det/mpasout.2024-05-02_13.00.00.nc' #1h forecast, init at 11utc, NO sfc continuous cycling
	STATIC='/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/251201-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L65_GFS'

	mkdir -p fig_horizontal_mpasjedi_diff
	python ./scripts-post/horizontal_mpasjedi_diff.py \
	    --filea  "$FILEA" \
	    --fileb  "$FILEB" \
	    --static "$STATIC"

fi




