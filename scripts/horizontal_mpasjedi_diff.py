#!/usr/bin/env python
from matplotlib.tri import Triangulation, TriAnalyzer
import warnings
import colormap
import numpy as np
import matplotlib.ticker as mticker
import matplotlib as mpl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import cartopy.crs as ccrs
# import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import matplotlib
matplotlib.use('agg')

warnings.filterwarnings('ignore')

# ------------ USER INPUT ----------------------------------------------------------
plot_var = "Increment"
decimal = 3            # number of decimals to round for text boxes

plot_box_width = 72.   # define size of plot domain (units: lat/lon degrees)
plot_box_height = 36.
cen_lat = 34.5
cen_lon = -97.5

# Determine extent for plot domain
half = plot_box_width / 2.
left = cen_lon - half
right = cen_lon + half
half = plot_box_height / 2.
bot = cen_lat - half
top = cen_lat + half


def plot_inc(var_inc, decimal=2):
    max_inc = np.around(np.max(var_inc), decimal)
    min_inc = np.around(np.min(var_inc), decimal)
    # decide the color contours based on the increment values
    clevmax = max((abs(max_inc), abs(min_inc)))
    inc = 0.1 * clevmax
    clevs = np.arange(-1.0 * clevmax, 1.0 * clevmax + inc, inc)
    cm = colormap.diff_colormap(clevs)
    return clevs, cm


def main():
    # janalysis = "/scratch1/BMC/wrfruc/jjhu/rundir/wrkflow-test/Btuning/2024050601_lbc/uv233/singleob_rh4rv0_avgheight_std14/mpasin.nc"
    # jbackgrnd = "/scratch1/BMC/wrfruc/jjhu/rundir/wrkflow-test/Btuning/2024050601_tuneB/bkg/mpasout.2024-05-06_01.00.00.nc"
    # jstatic = "/scratch1/BMC/wrfruc/jjhu/rundir/wrkflow-test/Btuning/2024050601_tuneB/invariant.nc"

    
    jfilea = '/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/250801-rrfs-workflow/OPSROOT/hrly_12km21/com/rrfs/v2.1.1/rrfs.20240527/00/ic/det/init.nc'
    jfileb = '/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/250801-rrfs-workflow/OPSROOT/hrly_12km21/com/rrfs/v2.1.1/rrfs.20240527/05/fcst/det/mpasout.2024-05-27_06.00.00.nc'
    jstatic = "/gpfs/f6/bil-pmp/scratch/Liaofan.Lin/250801-rrfs-workflow/rrfs-workflow/fix/meshes/conus12km.invariant.nc_L60_GFS"

    figdir = "./fig_horizontal_mpasjedi_diff/"

    # varible to plot
    variable = "SMOIS" # LANDMASK, ISLTYP, IVGTYP, SMOIS

    # target_lat, target_lon = 36.265, -95.145

    # Open NETCDF4 dataset for reading
    nc_a  = Dataset(jfilea, mode='r')
    nc_b = Dataset(jfileb, mode='r')
    f_latlon = Dataset(jstatic, "r")

    # read lat,lon information
    lats = np.array(f_latlon.variables['latCell'][:]) * 180.0 / np.pi  # Latitude of cells, rad
    lons0 = np.array(f_latlon.variables['lonCell'][:]) * 180.0 / np.pi  # Longitude of cells, rad
    lons = np.where(lons0 > 180.0, lons0 - 360.0, lons0)
    # z = f_latlon.variables['zgrid'][:]  # Geometric height of layer interfaces, m MSL

    # Grab variables
    if variable == "SMOIS":
        units = "m3/m3"
        jedi_a = nc_a.variables['smois'][0, :, :]
        jedi_b = nc_b.variables['smois'][0, :, :]     
        
    if variable == "LANDMASK":
        units = "-"
        jedi_a = nc_a.variables['landmask'][:] 
        jedi_b = nc_b.variables['landmask'][:]
               
    if variable == "ISLTYP":
        units = "-"
        jedi_a = nc_a.variables['isltyp'][:] 
        jedi_b = nc_b.variables['isltyp'][:]               
 
    if variable == "IVGTYP":
        units = "-"
        jedi_a = nc_a.variables['ivgtyp'][:] 
        jedi_b = nc_b.variables['ivgtyp'][:]                  

    jedi_inc_all = jedi_a - jedi_b  # the increment

    for lev in range(6, 9, 1):
        
        jedi_a_lev  = jedi_a[:, lev]
        jedi_b_lev  = jedi_b[:, lev]
        jedi_diff   = jedi_a_lev - jedi_b_lev

        # Print some final stats
        print(f"Level: {lev+1},  max: {np.around(np.max(jedi_diff), decimal)}, min: {np.around(np.min(jedi_diff), decimal)}")

        
        ## CREATE PLOT ##############################
        fig = plt.figure(figsize=(12.5, 10))
            
        # use scatter to plot the field
        markersize=0.6
        alpha=0.9
        
        # use a LambertConformal projection, choose trulat1 and trulat2 around the center of your domain to reduce distortion.
        ref_lon=-97.5  # reference or central lon
        ref_lat=36.0   # reference or centeal lat
        trulat1=36.0   # first standard parallel, no distortion
        trulat2=36.0   # second standar parallel, no distortion
        
        projection = ccrs.LambertConformal(central_longitude=ref_lon, central_latitude=ref_lat, standard_parallels=(trulat1, trulat2))
        ax = plt.axes(projection=projection)
        
        # add costlines, country borders, state/province borders
        ax.coastlines(resolution='50m')  # '110m', '50m', or '10m'
        ax.add_feature(cfeature.BORDERS, linewidth=0.7)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
           
        
        # Get the colormap setting
        clevs, cm = plot_inc(jedi_diff)   
        
        cmap1=plt.cm.get_cmap('bwr_r',20)
        cmaplist = [cmap1(i) for i in range(cmap1.N)]
        cmaplist[9] = [1,1,1,1]
        cmaplist[10] = [1,1,1,1]
        cmap1 = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap1.N)
        
        
        sc = ax.scatter(lons, lats, c=jedi_diff, cmap=cmap1, \
                        vmin=clevs[0], vmax=clevs[-1], \
                        transform=ccrs.PlateCarree(), s=markersize, alpha=alpha)
        
        
        # Colorbar
        cbar = plt.colorbar(sc, label='', orientation='horizontal', shrink=0.8, aspect=50, pad=0.01)    
        cbar.set_label(units, size=10)
        cbar.ax.tick_params(labelsize=10, rotation=30)
        cbar.set_ticks(clevs)
        
        
        # Add titles, text, and save the figure
        plt.suptitle(f"{variable} Diff at Level: {lev+1}", fontsize=20, y=0.95)
        #subtitle1_minmax = f"min: {np.around(np.min(jedi_inc), decimal)}\nmax: {np.around(np.max(jedi_inc), decimal)}"
        #m1.text(left * 0.99, bot * 1.01, f"{subtitle1_minmax}", fontsize=6, ha='left', va='bottom')
        plt.tight_layout()
        plt.savefig(f"{figdir}/diff_{variable}_z{lev}.png", dpi=250, bbox_inches='tight')
        plt.close()

        # CREATE PLOT ##############################
        fig = plt.figure(figsize=(12.5, 10))
            
        # use scatter to plot the field
        markersize=0.6
        alpha=0.9
        
        # use a LambertConformal projection, choose trulat1 and trulat2 around the center of your domain to reduce distortion.
        ref_lon=-97.5  # reference or central lon
        ref_lat=36.0   # reference or centeal lat
        trulat1=36.0   # first standard parallel, no distortion
        trulat2=36.0   # second standar parallel, no distortion
        
        projection = ccrs.LambertConformal(central_longitude=ref_lon, central_latitude=ref_lat, standard_parallels=(trulat1, trulat2))
        ax = plt.axes(projection=projection)
        
        # add costlines, country borders, state/province borders
        ax.coastlines(resolution='50m')  # '110m', '50m', or '10m'
        ax.add_feature(cfeature.BORDERS, linewidth=0.7)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')

        if variable == "SMOIS":
            clevs = np.linspace(0, 0.5, 21)
        else:
            clevs = np.linspace(jedi_a_lev.min(), jedi_a_lev.max(), 21)
        
        # colormap
        cmap1=plt.cm.get_cmap('jet',20)
           
        sc = ax.scatter(lons, lats, c=jedi_a_lev, cmap=cmap1, \
                        vmin=clevs[0], vmax=clevs[-1],\
                        transform=ccrs.PlateCarree(), s=markersize, alpha=alpha)
        
        # Colorbar
        cbar = plt.colorbar(sc, label='', orientation='horizontal', shrink=0.8, aspect=50, pad=0.01)    
        cbar.set_label(units, size=10)
        cbar.ax.tick_params(labelsize=10, rotation=30)
        cbar.set_ticks(clevs)
        
        # Add titles, text, and save the figure
        plt.suptitle(f"{variable} File A at Level: {lev+1}", fontsize=20, y=0.95)
        #subtitle1_minmax = f"min: {np.around(np.min(jedi_inc), decimal)}\nmax: {np.around(np.max(jedi_inc), decimal)}"
        #m1.text(left * 0.99, bot * 1.01, f"{subtitle1_minmax}", fontsize=6, ha='left', va='bottom')
        plt.tight_layout()
        plt.savefig(f"{figdir}/file_a_{variable}_z{lev}.png", dpi=250, bbox_inches='tight')
        plt.close()
        
                
        # CREATE PLOT ##############################
        fig = plt.figure(figsize=(12.5, 10))
            
        # use scatter to plot the field
        markersize=0.6
        alpha=0.9
        
        # use a LambertConformal projection, choose trulat1 and trulat2 around the center of your domain to reduce distortion.
        ref_lon=-97.5  # reference or central lon
        ref_lat=36.0   # reference or centeal lat
        trulat1=36.0   # first standard parallel, no distortion
        trulat2=36.0   # second standar parallel, no distortion
        
        projection = ccrs.LambertConformal(central_longitude=ref_lon, central_latitude=ref_lat, standard_parallels=(trulat1, trulat2))
        ax = plt.axes(projection=projection)
        
        # add costlines, country borders, state/province borders
        ax.coastlines(resolution='50m')  # '110m', '50m', or '10m'
        ax.add_feature(cfeature.BORDERS, linewidth=0.7)
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')
        
        if variable == "SMOIS":
            clevs = np.linspace(0, 0.5, 21)
        else:
            clevs = np.linspace(jedi_a_lev.min(), jedi_a_lev.max(), 21)
        
        # colormap
        cmap1=plt.cm.get_cmap('jet',20)
           
                      
        sc = ax.scatter(lons, lats, c=jedi_b_lev, cmap=cmap1, \
                        vmin=clevs[0], vmax=clevs[-1],\
                        transform=ccrs.PlateCarree(), s=markersize, alpha=alpha)
        
        # Colorbar
        cbar = plt.colorbar(sc, label='', orientation='horizontal', shrink=0.8, aspect=50, pad=0.01)    
        cbar.set_label(units, size=10)
        cbar.ax.tick_params(labelsize=10, rotation=30)
        cbar.set_ticks(clevs)
        
        # Add titles, text, and save the figure
        plt.suptitle(f"{variable} File B at Level: {lev+1}", fontsize=20, y=0.95)
        #subtitle1_minmax = f"min: {np.around(np.min(jedi_inc), decimal)}\nmax: {np.around(np.max(jedi_inc), decimal)}"
        #m1.text(left * 0.99, bot * 1.01, f"{subtitle1_minmax}", fontsize=6, ha='left', va='bottom')
        plt.tight_layout()
        plt.savefig(f"{figdir}/file_b_{variable}_z{lev}.png", dpi=250, bbox_inches='tight')
        plt.close()
        
        
                        
        
if __name__ == '__main__':
    main()
