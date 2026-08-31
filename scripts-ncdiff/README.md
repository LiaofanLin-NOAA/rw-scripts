# Check NetCDF Differences

This script checks and summarizes differences between two NetCDF files.

## Usage

### 1. Link the two NetCDF files

Create symbolic links to the two NetCDF files to be compared in the
`scripts-check-ncdiff` directory.

For example:

``` bash
ln -s /path/to/file1.nc file1.nc
ln -s /path/to/file2.nc file2.nc
```

### 2. Generate the NetCDF difference file

Load the NCO module and use `ncdiff` to generate a NetCDF file
containing the differences:

``` bash
module load nco
ncdiff file1.nc file2.nc ncdiff.nc
```

### 3. Load the Python environment

Load a Python environment containing the required packages. For example:

``` bash
module load pyDAmonitor
```

### 4. Run the check

Run `check_ncdiff.py` using the NetCDF difference file:

``` bash
python check_ncdiff.py ncdiff.nc
```
