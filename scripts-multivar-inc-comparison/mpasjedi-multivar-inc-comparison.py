#!/usr/bin/env python

import warnings
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from netCDF4 import Dataset


warnings.filterwarnings("ignore")


# =============================================================================
# USER SETTINGS
# =============================================================================

figdir = "./figures"
decimal = 3

markersize = 0.6
alpha = 0.9

# Lambert Conformal projection
ref_lon = -97.5
ref_lat = 36.0
trulat1 = 36.0
trulat2 = 36.0

# Number of color intervals
ncolors = 20


# =============================================================================
# HELPERS
# =============================================================================

def create_map_axis():
    """Create the Lambert Conformal map axis."""

    projection = ccrs.LambertConformal(
        central_longitude=ref_lon,
        central_latitude=ref_lat,
        standard_parallels=(trulat1, trulat2),
    )

    ax = plt.axes(projection=projection)

    ax.coastlines(resolution="50m")
    ax.add_feature(cfeature.BORDERS, linewidth=0.7)
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor="gray")

    return ax


def read_variable(nc_a, nc_b, variable):
    """Read the requested variable from File A and File B."""

    variable = variable.lower()

    if variable == "smois":
        units = "m3/m3"
        jedi_a = np.asarray(nc_a.variables["smois"][0, :, :])
        jedi_b = np.asarray(nc_b.variables["smois"][0, :, :])

    elif variable == "tslb":
        units = "K"
        jedi_a = np.asarray(nc_a.variables["tslb"][0, :, :])
        jedi_b = np.asarray(nc_b.variables["tslb"][0, :, :])

    elif variable == 'theta':
        units = "K"
        jedi_a = nc_a.variables["theta"][0, :, :].astype(np.float64)  # (Time, nCells, nVertLevels)
        jedi_b = nc_b.variables["theta"][0, :, :].astype(np.float64)
        pres_a = (nc_a.variables['pressure_p'][0, :, :] + nc_b['pressure_base'][0, :, :])/100.0
        pres_b = (nc_b.variables['pressure_p'][0, :, :] + nc_b['pressure_base'][0, :, :])/100.0
        dividend_a = (1000.0/pres_a)**(0.286)
        dividend_b = (1000.0/pres_b)**(0.286)
        jedi_a = jedi_a / dividend_a
        jedi_b = jedi_b / dividend_b 

    elif variable == 'qv':
        units = "g/kg"
        jedi_a = nc_a.variables['qv'][0, :, :] * 1000.0
        jedi_b = nc_b.variables['qv'][0, :, :] * 1000.0

    else:
        raise ValueError(
            f"Unsupported variable: {variable}. "
            "Currently supported variables are: smois, tslb"
        )

    return jedi_a, jedi_b, units


def get_shared_absolute_scale(variable, field_a, field_b):
    """
    Return ONE common color scale for File A and File B.

    File A and File B will use the exact same:
      * levels
      * colormap
      * normalization
    """

    variable = variable.lower()

    if variable == "smois":
        # Fixed soil-moisture range
        vmin = 0.0
        vmax = 0.5

    else:
        # Use the overall min/max from BOTH files
        vmin = min(
            float(np.nanmin(field_a)),
            float(np.nanmin(field_b)),
        )

        vmax = max(
            float(np.nanmax(field_a)),
            float(np.nanmax(field_b)),
        )

    levels = np.linspace(vmin, vmax, ncolors + 1)

    cmap = mpl.colormaps["jet"].resampled(ncolors)

    norm = mpl.colors.BoundaryNorm(
        boundaries=levels,
        ncolors=cmap.N,
        clip=True,
    )

    return levels, cmap, norm


def get_difference_scale(diff):
    """Create a symmetric difference color scale centered on zero."""

    max_abs = float(np.nanmax(np.abs(diff)))

    if max_abs == 0.0:
        return None, None, None

    levels = np.linspace(
        -max_abs,
        max_abs,
        ncolors + 1,
    )

    base_cmap = mpl.colormaps["bwr_r"].resampled(ncolors)
    colors = [base_cmap(i) for i in range(base_cmap.N)]

    # Make the two bins nearest zero white
    colors[9] = [1.0, 1.0, 1.0, 1.0]
    colors[10] = [1.0, 1.0, 1.0, 1.0]

    cmap = mpl.colors.ListedColormap(colors)

    norm = mpl.colors.BoundaryNorm(
        boundaries=levels,
        ncolors=cmap.N,
        clip=True,
    )

    return levels, cmap, norm


def plot_field(
    lons,
    lats,
    field,
    levels,
    cmap,
    norm,
    units,
    title,
    output_file,
):
    """Plot one field using an explicitly supplied color scale."""

    fig = plt.figure(figsize=(12.5, 10))
    ax = create_map_axis()

    sc = ax.scatter(
        lons,
        lats,
        c=field,
        cmap=cmap,
        norm=norm,
        transform=ccrs.PlateCarree(),
        s=markersize,
        alpha=alpha,
    )

    cbar = plt.colorbar(
        sc,
        orientation="horizontal",
        shrink=0.8,
        aspect=50,
        pad=0.01,
        boundaries=levels,
        ticks=levels,
    )

    cbar.set_label(units, size=10)
    cbar.ax.tick_params(labelsize=10, rotation=30)

    plt.suptitle(
        title,
        fontsize=16,
        y=0.98,
    )

    plt.tight_layout(
        rect=[0.0, 0.0, 1.0, 0.98]
    )

    plt.savefig(
        output_file,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main(filea, fileb, static, variable, level):

    variable = variable.lower()

    os.makedirs(figdir, exist_ok=True)

    print("")
    print("============================================================")
    print("MPAS-JEDI multivariable comparison")
    print("============================================================")
    print("FILE A   :", filea)
    print("FILE B   :", fileb)
    print("STATIC   :", static)
    print("VARIABLE :", variable)
    print("LEVEL    :", level)
    print("")

    # -------------------------------------------------------------------------
    # Read files
    # -------------------------------------------------------------------------

    with Dataset(filea, mode="r") as nc_a, \
         Dataset(fileb, mode="r") as nc_b, \
         Dataset(static, mode="r") as f_latlon:

        lats = (
            np.asarray(f_latlon.variables["latCell"][:])
            * 180.0
            / np.pi
        )

        lons0 = (
            np.asarray(f_latlon.variables["lonCell"][:])
            * 180.0
            / np.pi
        )

        lons = np.where(
            lons0 > 180.0,
            lons0 - 360.0,
            lons0,
        )
        
        landmask = np.asarray(
            f_latlon.variables["landmask"][:]
        )

        jedi_a, jedi_b, units = read_variable(
            nc_a,
            nc_b,
            variable,
        )

    # -------------------------------------------------------------------------
    # Convert 1-based level to 0-based Python index
    #
    # --level 6 means Python index 5.
    # -------------------------------------------------------------------------

    lev = level - 1

    if lev < 0 or lev >= jedi_a.shape[1]:
        raise ValueError(
            f"Requested level {level} is outside the available range "
            f"1-{jedi_a.shape[1]} for variable '{variable}'."
        )

    jedi_a_lev = np.asarray(jedi_a[:, lev])
    jedi_b_lev = np.asarray(jedi_b[:, lev])
    jedi_diff = jedi_a_lev - jedi_b_lev

    # -------------------------------------------------------------------------
    # Difference statistics over LAND only
    # landmask = 1: land
    # landmask = 0: water
    # -------------------------------------------------------------------------

    land_diff = jedi_diff[landmask == 1]

    land_count = land_diff.size
    land_min   = float(np.nanmin(land_diff))
    land_max   = float(np.nanmax(land_diff))
    land_mean  = float(np.nanmean(land_diff))
    land_std   = float(np.nanstd(land_diff))

    # -------------------------------------------------------------------------
    # Statistics
    # -------------------------------------------------------------------------

    a_min = float(np.nanmin(jedi_a_lev))
    a_max = float(np.nanmax(jedi_a_lev))

    b_min = float(np.nanmin(jedi_b_lev))
    b_max = float(np.nanmax(jedi_b_lev))

    diff_min = float(np.nanmin(jedi_diff))
    diff_max = float(np.nanmax(jedi_diff))

    print(
        f"File A level {level}: "
        f"min={a_min:.{decimal}f}, "
        f"max={a_max:.{decimal}f}"
    )

    print(
        f"File B level {level}: "
        f"min={b_min:.{decimal}f}, "
        f"max={b_max:.{decimal}f}"
    )

    print(
        f"Difference level {level}: "
        f"min={diff_min:.{decimal}f}, "
        f"max={diff_max:.{decimal}f}"
    )

    # =========================================================================
    # ONE SHARED COLOR SCALE FOR FILE A AND FILE B
    # =========================================================================

    absolute_levels, absolute_cmap, absolute_norm = \
        get_shared_absolute_scale(
            variable,
            jedi_a_lev,
            jedi_b_lev,
        )

    print("")
    print("Shared File A / File B color scale:")
    print(f"  min = {absolute_levels[0]:.{decimal}f}")
    print(f"  max = {absolute_levels[-1]:.{decimal}f}")
    print("")

    # =========================================================================
    # DIFFERENCE COLOR SCALE
    # =========================================================================

    diff_levels, diff_cmap, diff_norm = \
        get_difference_scale(jedi_diff)

    # =========================================================================
    # FILENAMES
    #
    # Example for TSLB level 6:
    #
    #   tslb_z6_fileA.png
    #   tslb_z6_fileB.png
    #   tslb_z6_diff.png
    # =========================================================================

    file_a_png = os.path.join(
        figdir,
        f"{variable}_z{level}_fileA.png",
    )

    file_b_png = os.path.join(
        figdir,
        f"{variable}_z{level}_fileB.png",
    )

    diff_png = os.path.join(
        figdir,
        f"{variable}_z{level}_diff.png",
    )

    # =========================================================================
    # FILE A
    # =========================================================================

    plot_field(
        lons=lons,
        lats=lats,
        field=jedi_a_lev,
        levels=absolute_levels,
        cmap=absolute_cmap,
        norm=absolute_norm,
        units=units,
        title=f"{variable.upper()} File A (Analysis) at Level: {level}; Unit: {units}\n",
        output_file=file_a_png,
    )

    # =========================================================================
    # FILE B
    #
    # Uses the EXACT SAME absolute_levels, absolute_cmap, and absolute_norm
    # as File A.
    # =========================================================================

    plot_field(
        lons=lons,
        lats=lats,
        field=jedi_b_lev,
        levels=absolute_levels,
        cmap=absolute_cmap,
        norm=absolute_norm,
        units=units,
        title=f"{variable.upper()} File B (Background) at Level: {level}; Unit: {units}\n",
        output_file=file_b_png,
    )

    # =========================================================================
    # DIFFERENCE
    # =========================================================================

    if diff_levels is None:
        print(
            "Warning: difference between File A and File B "
            "is zero everywhere."
        )
        print("No difference plot will be generated.")

    else:
        plot_field(
            lons=lons,
            lats=lats,
            field=jedi_diff,
            levels=diff_levels,
            cmap=diff_cmap,
            norm=diff_norm,
            units=units,
            title=(
                f"{variable.upper()} Diff (Analysis Increment) at Level: {level}; Unit: {units}\n"
                f"Over land: n={land_count}, "
                f"max={np.nanmax(land_diff):.{decimal}f}, "
                f"min={np.nanmin(land_diff):.{decimal}f}, "
                f"mean={land_mean:.{decimal}f}, "
                f"std={land_std:.{decimal}f}"
            ),
            output_file=diff_png,
        )

    print("")
    print("Generated figures:")
    print(f"  {file_a_png}")
    print(f"  {file_b_png}")

    if diff_levels is not None:
        print(f"  {diff_png}")

    print("")


# =============================================================================
# COMMAND LINE
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Compare one MPAS variable and vertical level "
            "between two MPAS files."
        )
    )

    parser.add_argument(
        "--filea",
        required=True,
        help="Path to File A",
    )

    parser.add_argument(
        "--fileb",
        required=True,
        help="Path to File B",
    )

    parser.add_argument(
        "--static",
        required=True,
        help="Path to static MPAS mesh / metadata file",
    )

    parser.add_argument(
        "--variable",
        required=True,
        help="Variable to process, e.g. smois or tslb",
    )

    parser.add_argument(
        "--level",
        required=True,
        type=int,
        help="Vertical level using 1-based numbering",
    )

    args = parser.parse_args()

    main(
        filea=args.filea,
        fileb=args.fileb,
        static=args.static,
        variable=args.variable,
        level=args.level,
    )
