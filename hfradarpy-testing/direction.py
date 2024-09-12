
# This file is derived from the hfradarpy github repository. It is for plotting
# the direction of the vector to or from the antenna (marked by different
# colors on the graph). I separated this file for easier access and to have
# a faster run time. 
# Look at their notebook for more information and comments.

# https://github.com/rucool/hfradarpy/blob/master/examples/notebooks/plot_radials_and_totals.ipynb


from hfradarpy.radials import Radial, concat
import glob
import os
import xarray as xr

site = 'MARA/'

# Path to radial directory
radial_dir = '../radial-data/processed/' + site

# Use glob to find radial files (*
#files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

#ds = concat(sorted(files), method='gridded', enhance=True)
#tds = ds.isel(time=-1)

radial_file = radial_dir + 'RDLm_MARA_2024_04_29_1600.ruv'

r = Radial(radial_file)

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

extent = [tds.lon.min()-dx, tds.lon.max()+dx, tds.lat.min()-dy, tds.lat.max()+dy]
# extent = [-83.0259, -78.9407, 22.8575, 26.6224]

#print(extent)

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

ax.set_extent(extent)
#ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.LAND)

import numpy as np

time = tds.time
lon = tds.coords['lon'].data
lat = tds.coords['lat'].data
u = tds['u'].data
v = tds['v'].data

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

from matplotlib import colors

"""
Motion displays the direction (towards or away) from radar
"""
title = 'Radial Map: Towards/Away from radar'
cmap= 'bwr'

velocity = tds.velocity
velocity_temp = velocity.where(velocity > 0, other=-1)  # Going away from radar
color_clipped = velocity_temp.where(velocity < 0, other=1).data  # Going towards radar
offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

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
