# Multivariable Increment Comparison

Scripts for comparing multiple variables between two MPAS-JEDI output files.

## Files

- `driver-multivar-inc.sh` - Driver script for setting up and running the comparison.
- `mpasjedi-multivar-inc-comparison.py` - Python script for comparing variables and generating plots.
- `colormap.py` - Colormap utilities used by the plotting script.

## Usage

Set the paths to the two MPAS output files and the static mesh file in `driver-multivar-inc.sh`:

```bash
FILEA="path/to/file_a.nc"
FILEB="path/to/file_b.nc"
STATIC="path/to/invariant.nc"
