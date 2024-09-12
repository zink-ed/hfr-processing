
# This file is derived from the hfradarpy github repository. It is for plotting
# if the vector passed or failed the quality tests (marked by different colors 
# on the graph). I separated this file to test if the quality control tests
# were successful in filtering out data.
# Look at their notebook for more information and comments.

# https://github.com/rucool/hfradarpy/blob/master/examples/notebooks/plot_radials_and_totals.ipynb


from hfradarpy.radials import Radial, qc_radial_file
import glob
import os
import xarray as xr

site = 'MARA/'

# Path to radial directory
radial_dir = '../radial-data/raw/' + site
#radial_dir = '../radial-data/processed/' + site

# Use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

r = Radial(files[10])

tds = r.to_xarray('gridded', enhance=True).squeeze()
tds = tds.squeeze()

receiver_location = [float(x) for x in tds.Origin.split('  ')]
receiver_location.reverse()
print(receiver_location)

# Import matplotlib and cartopy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

extent = []

dx = dy = 0.2  # Area around the point of interest.

extent = [tds.lon.min()+0.5, tds.lon.max() + 0.1, tds.lat.min(), tds.lat.max()-1]

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

ax.set_extent(extent)
#ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.LAND)

import numpy.ma as ma
import numpy as np

time = tds.time
lon = tds.coords['lon'].data
lat = tds.coords['lat'].data
u = tds['u'].data
v = tds['v'].data

u = ma.masked_invalid(u)
v = ma.masked_invalid(v)

from oceans.ocfis import uv2spdir, spdir2uv

angle, speed = uv2spdir(u, v)  # convert u/v to angle and speed

u, v = spdir2uv(  # convert angle and speed back to u/v, normalizing the arrow sizes
    np.ones_like(speed),
    angle,
    deg=True
)

from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm, Normalize

scale=None
headwidth=3
headlength=5
headaxislength=5
sub=1
velocity_min = -40
velocity_max = 40
cbar_step = 10
offset = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)

color_clipped = tds.vector_flag.where(tds.vector_flag == 0, other=-1).data  # PRIM == 1 where vectors pass qc

offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
title = "Radial Map: QC Pass/Fail"
cmap = colors.ListedColormap(['red', 'blue'])

# Plot title
plt.title(f'{title}\n{tds.time.data}')

qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# plot arrows over pcolor
h = ax.quiver(
    lon[::sub],
    lat[::sub],
    u[::sub],
    v[::sub],
    color_clipped,
    **qargs
)

plt.show()
