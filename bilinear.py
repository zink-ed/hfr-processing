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
r = Radial(files[80])

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

    for r in range(rows):

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

def bilinear(l, l2, r, c):
    b00 = l[r + 1][c - 1]
    b01 = l[r - 1][c - 1]
    b10 = l[r + 1][c + 1]
    b11 = l[r - 1][c + 1]
    
    b1 = (b00 + b10) / 2
    b2 = (b01 + b11) / 2
    
    result = (b1 + b2) / 2
    
    l[r][c] = result
    
    l2.append(result)
    

def interpolation():
    for r in range(1, rows - 1):
        for c in range(1, cols - 1):
            if interpolated_matrix[r][c] != 1:
                bilinear(u_matrix, u_interpolated, r, c)
                bilinear(v_matrix, v_interpolated, r, c)
                lon_interpolated.append(lon_matrix[r][c])
                lat_interpolated.append(lat_matrix[r][c])
                interpolated_matrix[r][c] = 3

interpolation()            
            
'''
boundary = [0, 0, 0, 0]


def bi_interpolate(dist, r, c):
    for d in range(1, dist):
        if interpolated_matrix[r - d][c] == 0:
            boundary[0] = 1
        if interpolated_matrix[r + d][c] == 0:
            boundary[0] = 1
        if interpolated_matrix[r][c - d] == 0:
            boundary[0] = 1
        if interpolated_matrix[r][c + d] == 0:
            boundary[0] = 1           

def interpolation(max_dist):
    for r in range(rows):
        for c in range(cols):
            if interpolated_matrix[r][c] == 0:
                bi_interpolate(dist, r, c)



def interpolate(rdist, cdist, r, c):
    
    dist = max(rdist, cdist)
    
    increment = (lon_matrix[r][c] - lon_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        lon_interpolated.append(lon_matrix[r][c] - d * increment)
        if rdist > cdist:
            lon_matrix[r - d][c] = lon_interpolated[-1]
            interpolated_matrix[r - d][c] = 2
        else:
            lon_matrix[r][c - d] = lon_interpolated[-1]
            interpolated_matrix[r][c - d] = 3
    
    increment = (lat_matrix[r][c] - lat_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        lat_interpolated.append(lat_matrix[r][c] - d * increment)
        if rdist > cdist:
            lat_matrix[r - d][c] = lat_interpolated[-1]
        else:
            lat_matrix[r][c - d] = lat_interpolated[-1]
    
    increment = (u_matrix[r][c] - u_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        u_interpolated.append(u_matrix[r][c] - d * increment)
        if rdist > cdist:
            u_matrix[r - d][c] = u_interpolated[-1]
        else:
            u_matrix[r][c - d] = u_interpolated[-1]
    
    increment = (v_matrix[r][c] - v_matrix[r - rdist][c - cdist]) / dist
    #print(increment)
    for d in range(1, dist):
        v_interpolated.append(v_matrix[r][c] - d * increment)
        if rdist > cdist:
            v_matrix[r - d][c] = v_interpolated[-1]
        else:
            v_matrix[r][c - d] = v_interpolated[-1]

def interpolation_bearing(max_dist):
    for r in range(rows):
        dist = max_dist + 1
        for c in range(cols):
            if interpolated_matrix[r][c] == 1:
                if dist <= max_dist and dist > 0:
                    interpolate(0, dist + 1, r, c)
                dist = 0
            else:
                dist = dist + 1
                
def interpolation_range(max_dist):
    for c in range(cols):
        dist = max_dist + 1
        for r in range(rows):
            if interpolated_matrix[r][c] == 1:
                if dist <= max_dist and dist > 0:
                    interpolate(dist + 1, 0, r, c)
                dist = 0
            else:
                dist = dist + 1     
'''
                
#interpolation_bearing(3)
#bearing_length = len(lon_interpolated)
#interpolation_range(3)                            

        
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

'''
colors = np.concatenate([['c'] * len(lon_original), ['tomato'] * bearing_length, ['orange'] * (len(lon_interpolated) - bearing_length)])

import matplotlib.lines as mlines
legend_lines = [mlines.Line2D([], [], color='c', lw=3, label='Original Data'), 
                 mlines.Line2D([], [], color='tomato', lw=3, label='Interpolated Bearing Data'),
                 mlines.Line2D([], [], color='orange', lw=3, label='Interpolated Range Data')]
'''

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
#ax.legend(handles=legend_lines)
plt.show()

