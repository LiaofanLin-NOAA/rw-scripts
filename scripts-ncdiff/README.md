# Check NetCDF Differences

This script checks and summarizes differences between two NetCDF files.

## Usage

### 1. Link the two NetCDF files

Create symbolic links to the two NetCDF files to be compared in the
`scripts-ncdiff` directory.

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

### 5. If the NetCDF files have differences

The script prints the full list of variables with nonzero differences,
including the minimum, maximum, and maximum absolute difference for each
variable.

For example:

``` text
Checking: diff-case04-mem-01-02.nc

acrunoff                            min=-1.477356e+01  max= 1.405677e+01  maxabs= 1.477356e+01
acsnom                              min=-8.255101e+00  max= 6.167744e+00  maxabs= 8.255101e+00
br                                  min=-2.023225e+00  max= 2.021233e+00  maxabs= 2.023225e+00
canwat                              min=-5.000269e-01  max= 5.002730e-01  maxabs= 5.002730e-01
cd                                  min=-2.840136e-02  max= 2.855286e-02  maxabs= 2.855286e-02
cda                                 min=-3.026308e-02  max= 3.046569e-02  maxabs= 3.046569e-02
ch                                  min=-1.603874e-02  max= 1.229774e-02  maxabs= 1.603874e-02
chklowq                             min=-1.000000e+00  max= 1.000000e+00  maxabs= 1.000000e+00
chs                                 min=-3.728140e-02  max= 4.405760e-02  maxabs= 4.405760e-02
chs2                                min=-4.667264e-02  max= 5.470352e-02  maxabs= 5.470352e-02
ck                                  min=-1.500636e-02  max= 1.171858e-02  maxabs= 1.500636e-02
cka                                 min=-1.607547e-02  max= 1.248822e-02  maxabs= 1.607547e-02
cldfrac                             min=-1.000000e+00  max= 1.000000e+00  maxabs= 1.000000e+00
cldfrac_bl                          min=-1.000000e+00  max= 1.000000e+00  maxabs= 1.000000e+00
cov                                 min=-1.119780e-02  max= 6.707015e-03  maxabs= 1.119780e-02
cpm                                 min=-6.626038e+00  max= 6.769348e+00  maxabs= 6.769348e+00
cqs                                 min=-3.732865e-02  max= 4.415457e-02  maxabs= 4.415457e-02
cqs2                                min=-4.700988e-02  max= 5.516101e-02  maxabs= 5.516101e-02
cuprec                              min=-1.806948e-02  max= 8.015150e-03  maxabs= 1.806948e-02
dew                                 min=-2.036965e-05  max= 1.750626e-05  maxabs= 2.036965e-05
dtaux3d                             min=-1.017919e-02  max= 1.249226e-02  maxabs= 1.249226e-02
dtauy3d                             min=-1.500356e-02  max= 1.622767e-02  maxabs= 1.622767e-02
dusfcg                              min=-1.433282e+00  max= 1.278419e+00  maxabs= 1.433282e+00
dvsfcg                              min=-2.827458e+00  max= 1.015671e+00  maxabs= 2.827458e+00
el_pbl                              min=-3.763976e+02  max= 3.788775e+02  maxabs= 3.788775e+02
exner                               min=-1.846254e-03  max= 1.981318e-03  maxabs= 1.981318e-03
flhc                                min=-4.089281e+01  max= 5.283684e+01  maxabs= 5.283684e+01
flqc                                min=-4.089281e+01  max= 5.283684e+01  maxabs= 5.283684e+01
glw                                 min=-8.333960e+01  max= 7.885745e+01  maxabs= 8.333960e+01
graupelnc                           min=-3.106751e-01  max= 5.446755e+00  maxabs= 5.446755e+00
graupelncv                          min=-5.542914e-03  max= 4.887026e-02  maxabs= 4.887026e-02
grdflx                              min=-1.402903e+02  max= 2.580338e+02  maxabs= 2.580338e+02
gsw                                 min=-4.005192e+02  max= 2.790573e+02  maxabs= 4.005192e+02
gz1oz0                              min=-2.742935e+00  max= 1.138336e+00  maxabs= 2.742935e+00
hfx                                 min=-1.652376e+02  max= 2.239264e+02  maxabs= 2.239264e+02
hpbl                                min=-3.776289e+03  max= 4.044631e+03  maxabs= 4.044631e+03
kpbl                                min=-2.200000e+01  max= 2.400000e+01  maxabs= 2.400000e+01
lh                                  min=-5.156570e+02  max= 6.080202e+02  maxabs= 6.080202e+02
mavail                              min=-7.343442e-01  max= 9.999900e-01  maxabs= 9.999900e-01
mol                                 min=-2.963149e-01  max= 2.509315e-01  maxabs= 2.963149e-01
o3vmr                               min=-4.107960e-08  max= 3.357945e-08  maxabs= 4.107960e-08
precipfr                            min=-2.165140e-02  max= 8.882255e-02  maxabs= 8.882255e-02
precipw                             min=-2.487036e+01  max= 2.516358e+01  maxabs= 2.516358e+01
pressure_p                          min=-3.815791e+02  max= 8.659487e+02  maxabs= 8.659487e+02
psih                                min=-2.087202e+01  max= 1.993115e+01  maxabs= 2.087202e+01
psim                                min=-2.319283e+01  max= 2.265427e+01  maxabs= 2.319283e+01
q2                                  min=-7.897751e-03  max= 8.249787e-03  maxabs= 8.249787e-03
qc                                  min=-1.759948e-03  max= 1.574075e-03  maxabs= 1.759948e-03
qc_bl                               min=-3.327213e-04  max= 3.673816e-04  maxabs= 3.673816e-04
qcg                                 min=-7.937793e-04  max= 8.791033e-04  maxabs= 8.791033e-04
qfx                                 min=-2.062628e-04  max= 2.432081e-04  maxabs= 2.432081e-04
qg                                  min=-4.088837e-03  max= 3.444747e-03  maxabs= 4.088837e-03
qgh                                 min=-3.794441e-02  max= 1.497138e-02  maxabs= 3.794441e-02
qi                                  min=-1.164552e-03  max= 8.961902e-04  maxabs= 1.164552e-03
qi_bl                               min=-3.089454e-05  max= 3.893521e-05  maxabs= 3.893521e-05
qke                                 min=-2.751827e+01  max= 1.798555e+01  maxabs= 2.751827e+01
qke_adv                             min=-2.751827e+01  max= 1.798555e+01  maxabs= 2.751827e+01
qr                                  min=-2.785420e-03  max= 2.516411e-03  maxabs= 2.785420e-03
qs                                  min=-8.445131e-03  max= 6.763687e-03  maxabs= 8.445131e-03
qsfc                                min=-1.904667e-02  max= 1.387813e-02  maxabs= 1.904667e-02
qsg                                 min=-3.029517e-02  max= 1.336701e-02  maxabs= 3.029517e-02
qsq                                 min=-1.653985e-06  max= 3.641675e-06  maxabs= 3.641675e-06
qv                                  min=-1.142582e-02  max= 9.879273e-03  maxabs= 1.142582e-02
qvg                                 min=-2.081909e-02  max= 1.418534e-02  maxabs= 2.081909e-02
rainc                               min=-3.585848e+01  max= 2.632437e+01  maxabs= 3.585848e+01
raincv                              min=-1.084169e+00  max= 4.809090e-01  maxabs= 1.084169e+00
rainncv                             min=-1.144905e+00  max= 8.491065e-01  maxabs= 1.144905e+00
re_cloud                            min=-4.736840e-05  max= 4.736840e-05  maxabs= 4.736840e-05
re_ice                              min=-1.125000e-04  max= 1.125000e-04  maxabs= 1.125000e-04
re_snow                             min=-9.990000e-04  max= 9.990000e-04  maxabs= 9.990000e-04
refl10cm                            min=-8.631071e+01  max= 8.657687e+01  maxabs= 8.657687e+01
refl10cm_max                        min=-6.875011e+01  max= 8.057731e+01  maxabs= 8.057731e+01
relhum                              min=-1.128263e+02  max= 1.480345e+02  maxabs= 1.480345e+02
rho                                 min=-3.374290e-02  max= 5.536079e-02  maxabs= 5.536079e-02
rho_p                               min=-2.411079e-02  max= 5.872995e-02  maxabs= 5.872995e-02
rho_zz                              min=-2.411079e-02  max= 5.873001e-02  maxabs= 5.873001e-02
rhosnf                              min=-1.500000e+03  max= 1.500000e+03  maxabs= 1.500000e+03
rmol                                min=-2.235322e+00  max= 2.271960e+00  maxabs= 2.271960e+00
rt_diabatic_tend                    min=-4.263509e-02  max= 2.693125e-02  maxabs= 4.263509e-02
rtheta_p                            min=-9.789734e-01  max= 1.487091e+00  maxabs= 1.487091e+00
ru                                  min=-3.777900e+01  max= 4.110979e+01  maxabs= 4.110979e+01
ru_p                                min=-3.943282e-01  max= 5.078717e-01  maxabs= 5.078717e-01
rubldiff                            min=-1.017919e-02  max= 1.249226e-02  maxabs= 1.249226e-02
rvbldiff                            min=-1.500356e-02  max= 1.622767e-02  maxabs= 1.622767e-02
rw                                  min=-4.922846e+00  max= 3.370748e+00  maxabs= 4.922846e+00
rw_p                                min=-3.142642e-01  max= 3.930286e-01  maxabs= 3.930286e-01
sfc_albedo                          min=-5.432164e-02  max= 2.909322e-01  maxabs= 2.909322e-01
sfc_emiss                           min=-1.325023e-02  max= 7.297993e-03  maxabs= 1.325023e-02
sfcrunoff                           min=-1.477356e+01  max= 1.405677e+01  maxabs= 1.477356e+01
sh2o                                min=-2.870930e-01  max= 4.140000e-01  maxabs= 4.140000e-01
skintemp                            min=-1.602820e+01  max= 8.470642e+00  maxabs= 1.602820e+01
smfr3d                              min=-4.072585e-01  max= 4.599999e-01  maxabs= 4.599999e-01
smois                               min=-2.870815e-01  max= 4.140000e-01  maxabs= 4.140000e-01
smstav                              min=-1.271511e+02  max= 1.279615e+02  maxabs= 1.279615e+02
snow                                min=-1.273245e+01  max= 7.111123e+00  maxabs= 1.273245e+01
snowc                               min=-2.297817e-01  max= 1.472841e-01  maxabs= 2.297817e-01
snowfallac                          min=-3.826805e-03  max= 1.079577e-02  maxabs= 1.079577e-02
snowh                               min=-6.286321e-02  max= 3.402592e-02  maxabs= 6.286321e-02
snownc                              min=-3.693950e-01  max= 1.958695e-01  maxabs= 3.693950e-01
snowncv                             min=-1.973384e-02  max= 3.424690e-03  maxabs= 1.973384e-02
soilt1                              min=-4.476361e+01  max= 2.343826e+01  maxabs= 4.476361e+01
sr                                  min=-1.000000e+00  max= 1.000000e+00  maxabs= 1.000000e+00
sst                                 min=-1.771423e+01  max= 7.159271e+00  maxabs= 1.771423e+01
sstsk                               min=-1.771423e+01  max= 7.159271e+00  maxabs= 1.771423e+01
surface_pressure                    min=-3.822578e+02  max= 8.264609e+02  maxabs= 8.264609e+02
t2m                                 min=-1.604135e+01  max= 6.074951e+00  maxabs= 1.604135e+01
tend_sfc_pressure                   min=-5.022135e-01  max= 7.096354e-01  maxabs= 7.096354e-01
th2m                                min=-1.653848e+01  max= 6.513519e+00  maxabs= 1.653848e+01
theta                               min=-1.933273e+01  max= 3.167688e+01  maxabs= 3.167688e+01
theta_m                             min=-1.812650e+01  max= 3.018259e+01  maxabs= 3.018259e+01
tlag                                min=-8.548218e+00  max= 4.771484e+00  maxabs= 8.548218e+00
tmn                                 min=-8.548218e+00  max= 4.771484e+00  maxabs= 8.548218e+00
tslb                                min=-1.602820e+01  max= 8.470642e+00  maxabs= 1.602820e+01
tsnav                               min=-2.048694e+01  max= 1.713898e+01  maxabs= 2.048694e+01
tsq                                 min=-2.999964e+00  max= 2.999872e+00  maxabs= 2.999964e+00
tyear_mean                          min=-8.548218e+00  max= 4.771484e+00  maxabs= 8.548218e+00
u                                   min=-4.741153e+01  max= 3.945668e+01  maxabs= 4.741153e+01
u10                                 min=-1.251214e+01  max= 2.550872e+01  maxabs= 2.550872e+01
uReconstructMeridional              min=-3.070368e+01  max= 3.085041e+01  maxabs= 3.085041e+01
uReconstructZonal                   min=-3.831651e+01  max= 3.844680e+01  maxabs= 3.844680e+01
ust                                 min=-6.206410e-01  max= 1.025834e+00  maxabs= 1.025834e+00
ustm                                min=-6.206410e-01  max= 1.025834e+00  maxabs= 1.025834e+00
v10                                 min=-2.055692e+01  max= 9.459208e+00  maxabs= 2.055692e+01
w                                   min=-9.113595e+00  max= 9.215187e+00  maxabs= 9.215187e+00
wspd                                min=-6.569079e+00  max= 1.299213e+01  maxabs= 1.299213e+01
z0                                  min=-8.314536e-02  max= 8.477058e-02  maxabs= 8.477058e-02
znt                                 min=-8.314536e-02  max= 8.477058e-02  maxabs= 8.477058e-02
zol                                 min=-2.005431e+01  max= 2.039036e+01  maxabs= 2.039036e+01
...
```

The actual run will include the full printout for all variables with
nonzero differences.

### 6. If the NetCDF files have no differences

If no differences are found, no variables are printed after the file
name.

For example:

``` text
Checking: diff-mem01-case-03-04.nc
```

This indicates that no differences were found between the two NetCDF
files.
