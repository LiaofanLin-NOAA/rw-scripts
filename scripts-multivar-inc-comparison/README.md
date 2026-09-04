# Multivariable Increment Comparison

Scripts for comparing variables between two MPAS-JEDI output files and generating horizontal plots.

## Files

* `driver-multivar-inc.sh` - Slurm driver script.
* `mpasjedi-multivar-inc-comparison.py` - Comparison and plotting script.
* `colormap.py` - Colormap utilities.

## Usage

Set the input files and variables/levels in `driver-multivar-inc.sh`.

```bash
# NOTE:
#   Atmospheric variables: level 1 = bottom model level/layer.
#   Land-surface variables: level 1 = top soil level/layer.
VAR_LEVELS=(
    "smois:1"
    "tslb:1"
    "theta:1"
    "qv:1"
)
```

Submit with:

```bash
sbatch driver-multivar-inc.sh
```

Figures are saved in `./figures/` for File A, File B, and File A minus File B.
