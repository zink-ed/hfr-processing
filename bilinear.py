# import libraries
from hfradarpy.radials import Radial
import glob
import os
import xarray as xr

# read in the data of a radial file
radial_dir = './radials_clean/MARA/'

# save directory
save_dir = './post_int_images'

# use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

# select file
r = Radial(files[0])

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

lon_matrix = np.zeros((rows, cols))
lat_matrix = np.zeros((rows, cols))
u_matrix = np.zeros((rows, cols))
v_matrix = np.zeros((rows, cols))
interpolated_matrix = np.zeros((rows, cols))

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
    interpolated_matrix[r][b] = 1
    
    
# getting interpolated data

u_interpolated = []
v_interpolated = []

# interpolated lon and lat

lon_interpolated = []
lat_interpolated = []

from scipy.interpolate import interp1d

def split():
    
    
    i = 0
    start = 0

    for r in range(rows - 10):

        temp_lon = []
        temp_lat = []
        temp_bear = []

        while ranges[i] == ranges[start]:
            temp_lon.append(lon_original[i])
            temp_lat.append(lat_original[i])
            temp_bear.append(bearings[i])
            i = i + 1
            
            if i == len(ranges):
                break
            
        new_angles = np.arange(temp_bear[0], temp_bear[-1], 2)
        
        offset = int((temp_bear[0] - minc) / 2)
        
        if_lon = interp1d(np.array(temp_bear), np.array(temp_lon), kind='quadratic')
        i_lon = if_lon(new_angles)
        
        if_lat = interp1d(np.array(temp_bear), np.array(temp_lat), kind='quadratic')
        i_lat = if_lat(new_angles)
        
        c = offset
        
        for n, t in zip(i_lon, i_lat):
            
            if interpolated_matrix[r][c] != 1:
                lon_matrix[r][c] = n
                lat_matrix[r][c] = t
                #lon_interpolated.append(n)
                #lat_interpolated.append(t)
                interpolated_matrix[r][c] = 2
            c = c + 1
       
        start = i

split()

# POSSIBLE GRIDS

# 1: 1, 1, 1, 1

# 2: 1, 1, 1, 2
#    1, 1, 2, 1
#    1, 2, 1, 1
#    2, 1, 1, 1

#    2, 2, 1, 1
#    2, 1, 2, 1
#    2, 1, 1, 2
#    1, 2, 2, 1
#    1, 2, 1, 2
#    1, 1, 2, 2

#    2, 2, 2, 1
#    1, 2, 2, 2
#    2, 1, 2, 2
#    2, 2, 1, 2

#    2, 2, 2, 2

def check_grid(r, c, a):
    
    if (interpolated_matrix[r + a[0]][c - a[2]] != 1 or 
        interpolated_matrix[r - a[1]][c - a[2]] != 1 or
        interpolated_matrix[r + a[0]][c + a[3]] != 1 or
        interpolated_matrix[r - a[1]][c + a[3]] != 1):
            return False
    
    return True

def get_grid(r, c, a):

    if check_grid(r, c, a):
        return True

    for e in a:
        e = 2
        if check_grid(r, c, a):
            return True
        e = 1
        
    for i in range(4):
        a[i] = 2
        for j in range(i + 1, 4):
            a[j] = 2
            if check_grid(r, c, a):
                return True
            a[j] = 1
        a[i] = 1     
    
    return False
        


# array: [dr1, dr2, dc1, dc2]

def bilinear(m, l, r, c, a):

    b00 = m[r + a[0]][c - a[2]]
    b01 = m[r - a[1]][c - a[2]]
    b10 = m[r + a[0]][c + a[3]]
    b11 = m[r - a[1]][c + a[3]]
    
    b1 = (b00 * a[0] + b01 * a[1]) / (a[0] + a[1])
    b2 = (b10 * a[0] + b11 * a[1]) / (a[0] + a[1])
    
    result = (b1 * a[2] + b2 * a[3]) / (a[2] + a[3])
    
    m[r][c] = result
    l.append(result)


def interpolation(d):

    for r in range(d, rows - d):
    
        for c in range(d, cols - d):
            
            a = [1, 1, 1, 1]
            
            # if value has been interpolated
            if interpolated_matrix[r][c] == 2:
            
            # run through possible grids for bilinear
                if get_grid(r, c, a):
                    bilinear(u_matrix, u_interpolated, r, c, a)
                    bilinear(v_matrix, v_interpolated, r, c, a)
                    lon_interpolated.append(lon_matrix[r][c])
                    lat_interpolated.append(lat_matrix[r][c])
                    interpolated_matrix[r][c] = 3

#interpolation(1)        
interpolation(2)    

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

# plot details
title = 'Radial Map: Original and Interpolated Data'
plt.title(f'{title}\n')

offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# turning into numpy arrays (combine both)
x = np.concatenate((lon_original, lon_interpolated)) 
y = np.concatenate((lat_original, lat_interpolated)) 
u = np.concatenate((u_original, u_interpolated))
v = np.concatenate((v_original, v_interpolated)) 

qargs = dict(scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# add color

colors = np.concatenate([['c'] * len(lon_original), ['tomato'] * (len(lon_interpolated))])



import matplotlib.lines as mlines
legend_lines = [mlines.Line2D([], [], color='c', lw=3, label='Original Data'), 
                 mlines.Line2D([], [], color='tomato', lw=3, label='Interpolated Data')]


# plot arrows
q1 = ax.quiver(
    x,
    y,
    u,
    v,
    color=colors,
    **qargs
)



# include legend
#plt.scatter(x, y)
ax.legend(handles=legend_lines)
plt.show()

