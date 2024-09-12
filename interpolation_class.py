
# Created Aug 2024 by Cathleen Qiao

# This file was created to deal with interpolation issues. I had an idea to 
# create a class that I could use for interpolating the radial data instead 
# of creating so many matrices. 

# This file is used for testing.


from hfradarpy.radials import Radial
import xarray as xr

# read in the data of a radial file
radial_dir = './radials_clean/'
radial_file = 'RDLm_MARA_2024_04_29_1600.ruv'

r = Radial(radial_dir + radial_file)

# xarray
#tds = r.to_xarray('gridded', enhance=True).squeeze()
#tds = tds.squeeze()

# dataframe
df = r.data

# separate the data to have specific variables
antenna_lon = -80.9832833
antenna_lat = 24.7401333

lon = df['LOND'].to_numpy()
lat = df['LATD'].to_numpy()

rg = df['SPRC'].to_numpy()

uv = df['VELU'].to_numpy()
vv = df['VELV'].to_numpy()

ranges = df['RNGE'].to_numpy()
bearings = df['BEAR'].to_numpy()

vel = df['VELO'].to_numpy()

# creating class

class Data:
    def __init__(self, lon, lat, u, v, i = True):
        self.lon = lon
        self.lat = lat
        self.u = u
        self.v = v
        self.i = i

# create matrix

import numpy as np

# getting dimensions
maxr = np.max(rg)
minr = np.min(rg)
maxc = np.max(bearings)
minc = np.min(bearings)

rows = maxr - minr + 1
cols = int((maxc - minc) / 2) + 1

default = Data(0, 0, 0, 0)

#matrix = [[default for _ in range(cols)] for _ in range(rows)]
#matrix = np.zeros((maxr - minr + 1, int((maxc - minc) / 2) + 1))

lonm = np.zeros((rows, cols))
latm = np.zeros((rows, cols))
um = np.zeros((rows, cols))
vm = np.zeros((rows, cols))
im = np.full((rows, cols), True, dtype=bool)

# putting in radial velocity into matrix

for r, b, n, t, u, v in zip(rg, bearings, lon, lat, uv, vv):
    # r for the row
    r = r - minr
    # b for the column
    b = int((b - minc) / 2)
    # setting matrix element
    #matrix[r][b] = Data(n, t, u, v, False)
    lonm[r][b] = n
    latm[r][b] = t
    um[r][b] = u
    vm[r][b] = v
    im[r][b] = False
    
def get_arrays(matrix):
    n = []
    t = []
    u = []
    v = []
    i = []
    for r in matrix:
        for e in r:
            n.append(e.lon)
            t.append(e.lat)
            u.append(e.u)
            v.append(e.v)
            i.append(e.i)
    return np.array(n), np.array(t), np.array(u), np.array(v), i
    
# plotting from the matrix (to monitor progress)
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

extent = [-82.8, -78.3, 22.5, 26.1]

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

ax.set_extent(extent)
#ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.LAND)

import numpy as np

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

title = 'Radial Map: Towards/Away from Radar'
cmap= 'bwr'

#velocity = tds.velocity
#velocity_temp = velocity.where(velocity > 0, other=-1)  # Going away from radar
#color_clipped = velocity_temp.where(velocity < 0, other=1).data  # Going towards radar
offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# Plot title
plt.title(f'{title}\n')

qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

n, t, u, v, i = get_arrays(matrix)

#colors = ['blue' if it else 'green' for it in i] 

from matplotlib.colors import to_rgba

colors = np.where(im, to_rgba('green'), to_rgba('blue'))

# plot arrows over pcolor
h = ax.quiver(
    lonm,
    latm,
    um,
    vm,
    color=colors,
    **qargs
)

plt.show()

# create a new flag that shows interpolated data (use this to plot later)


# bearing

# line of best fit: quadratic?


# radial

# line of best fit: linear
