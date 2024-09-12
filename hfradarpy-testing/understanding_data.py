
# This file is for understanding the different types of data provided 
# in the radial files. I plotted the LOND and LATD and compared
# it with the XDST and YDST to see if there were any differences.

# Look at the hfradarpy repository for more information on the radial class:
# https://github.com/rucool/hfradarpy/blob/master/hfradarpy/radials.py


# import libraries
from hfradarpy.radials import Radial
import glob
import os
import xarray as xr

# read in the data of a radial file
radial_dir = '../radial-data/processed/MARA/'


# use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

r = Radial(files[10])

# dataframe
df = r.data

# separate the data to have specific variables (put into numpy arrays)
#antenna_lon = -80.9832833
#antenna_lat = 24.7401333

lon_original = df['LOND'].to_numpy()
lat_original = df['LATD'].to_numpy()

rg = df['SPRC'].to_numpy()

u_original = df['VELU'].to_numpy()
v_original = df['VELV'].to_numpy()

ranges = df['RNGE'].to_numpy()
bearings = df['BEAR'].to_numpy()

x = df['XDST'].to_numpy()
y = df['YDST'].to_numpy()

vel = df['VELO'].to_numpy()

# plot

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Intialize an empty subplot using cartopy
#fig, ax = plt.subplots(
#    figsize=(11, 8),
#    subplot_kw=dict(projection=ccrs.Mercator())
#)

# comparing lon and lat with x and y
fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 8))

# setting the extent
extent = [-82.8, -78.3, 22.5, 26.1]
#ax.set_extent(extent)

colors = ('b', 'r')

#ax.append(plt.axes())
ax1.scatter(lon_original, lat_original, color='blue')
ax1.autoscale()

ax2.scatter(x, y, color='red', alpha=0.7)
ax2.autoscale()
plt.show()
