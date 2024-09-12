# import libraries
from hfradarpy.radials import Radial
import xarray as xr

# read in the data of a radial file
radial_dir = './radials_clean/'
radial_file = 'RDLm_MARA_2024_04_29_1600.ruv'

r = Radial(radial_dir + radial_file)

# dataframe
df = r.data

# separate the data to have specific variables (put into numpy arrays)
antenna_lon = -80.9832833
antenna_lat = 24.7401333

lon_original = df['LOND'].to_numpy()
lat_original = df['LATD'].to_numpy()

rg = df['SPRC'].to_numpy()

u_original = df['VELU'].to_numpy()
v_original = df['VELV'].to_numpy()

ranges = df['RNGE'].to_numpy()
bearings = df['BEAR'].to_numpy()

vel = df['VELO'].to_numpy()

# create matrix

import numpy as np

# getting dimensions
maxr = np.max(rg)
minr = np.min(rg)
maxc = np.max(bearings)
minc = np.min(bearings)

rows = maxr - minr + 1
cols = int((maxc - minc) / 2) + 1

#matrix = np.zeros((maxr - minr + 1, int((maxc - minc) / 2) + 1))

lon_matrix = np.zeros((rows, cols))
lat_matrix = np.zeros((rows, cols))
u_matrix = np.zeros((rows, cols))
v_matrix = np.zeros((rows, cols))
interpolated_matrix = np.full((rows, cols), True, dtype=bool)

# putting in radial velocity into matrix

for r, b, n, t, u, v in zip(rg, bearings, lon_original, lat_original, u_original, v_original):
    # r for the row
    r = r - minr
    # b for the column
    b = int((b - minc) / 2)
    
    # setting matrix element
    lon_matrix[r][b] = n
    lat_matrix[r][b] = t
    u_matrix[r][b] = u
    v_matrix[r][b] = v
    interpolated_matrix[r][b] = False
    
    
# getting interpolated data

lon_interpolated = []
lat_interpolated = []
u_interpolated = []
v_interpolated = []

#g = int(input("Range: "))

def interpolate(rdist, cdist, r, c):
    
    dist = max(rdist, cdist)
    
    increment = (lon_matrix[r][c] - lon_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        lon_interpolated.append(d * increment + lon_matrix[r - rdist][c - cdist])
    
    increment = (lat_matrix[r][c] - lat_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        lat_interpolated.append(d * increment + lat_matrix[r - rdist][c - cdist])
    
    increment = (u_matrix[r][c] - u_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        u_interpolated.append(d * increment + u_matrix[r - rdist][c - cdist])
    
    increment = (v_matrix[r][c] - v_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        v_interpolated.append(d * increment + v_matrix[r - rdist][c - cdist])

def interpolation_bearing(max_dist):
    for r in range(rows):
        dist = max_dist + 1
        for c in range(cols):
            if interpolated_matrix[r][c] == False:
                if dist <= max_dist and dist > 0:
                    interpolate(0, dist + 1, r, c)
                dist = 0
            else:
                dist = dist + 1
                
def interpolation_range(max_dist):
    for c in range(cols):
        dist = max_dist + 1
        for r in range(rows):
            if interpolated_matrix[r][c] == False:
                if dist <= max_dist and dist > 0:
                    interpolate(dist + 1, 0, r, c)
                dist = 0
            else:
                dist = dist + 1
    
interpolation_bearing(3)
interpolation_range(3)

# saving interpolated data to csv
"""
import pandas as pd

interpolated_data = pd.DataFrame({'lon' : lon_interpolated,
                                  'lat' : lat_interpolated, 
                                  'u' : u_interpolated,
                                  'v' : v_interpolated })
                                  
interpolated_data.to_csv('interpolated.csv', sep=' ', index=False)     
"""                             
        
# plotting from the matrix (to monitor progress)
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

# setting the extent
extent = [-82.8, -78.3, 22.5, 26.1]
ax.set_extent(extent)

# import plotting stuff
from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm, Normalize

# quiver configurations
scale=None
headwidth=3
headlength=5
headaxislength=5
#sub=1
velocity_min = -40
velocity_max = 40
#cbar_step = 10
#offset = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)

# attempting to get scale
"""
def get_scale(n, t, c, d):
    q = plt.quiver(lon_matrix, lat_matrix, u_matrix, v_matrix)
    plt.close()
    return q.scale
"""

# plot details
title = 'Radial Map: Original and Interpolated Data'
plt.title(f'{title}\n')

offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# attempting to normalize data
"""
n = lon_matrix / np.linalg.norm(lon_matrix)
t = lat_matrix / np.linalg.norm(lat_matrix)
c = u_matrix / np.linalg.norm(u_matrix)
d = v_matrix / np.linalg.norm(v_matrix)
"""

# turning into numpy arrays
x = np.array(lon_interpolated) 
y = np.array(lat_interpolated) 
u = np.array(u_interpolated)
v = np.array(v_interpolated) 

qargs = dict(scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# plot arrows over pcolor
q1 = ax.quiver(
    lon_matrix,
    lat_matrix,
    u_matrix,
    v_matrix,
    color='blue',
    label='Original Data',
    **qargs
)

q2 = ax.quiver(
    x,
    y,
    u,
    v,
    color='green',
    label='Interpolated Data',
    **qargs
)

ax.legend()
plt.show()

