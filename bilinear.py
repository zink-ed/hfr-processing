# import libraries
from hfradarpy.radials import Radial
import glob
import os
import xarray as xr

# read in the data of a radial file
radial_dir = './radials_clean/MARA/'
#radial_dir = './radials_clean/KEYWEST/'

# save directory
save_dir = './post_int_images'

# use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

# select file
r = Radial(files[70])

# dataframe
df = r.data

# separate the data to have specific variables (put into numpy arrays)
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
    interpolated_matrix[r][b] = 1   # 1 means original
       
# getting interpolated u and v data
u_lin = []
v_lin = []
u_bilin = []
v_bilin = []

# interpolated lon and lat
lon_lin = []
lat_lin = []

lon_bilin = []
lat_bilin = []

# dictionaries
temp_lon = {}
temp_lat = {}
temp_bear = {}

# separate list into dict by range
def create_dict(list1, dict1):
    for r, l in zip(ranges, list1):
        if r in dict1:
            dict1[r].append(l)
        else:
            dict1[r] = [l]
        
create_dict(lon_original, temp_lon)
create_dict(lat_original, temp_lat)
create_dict(bearings, temp_bear)

#print(temp_lon)
#print(temp_lat)
#print(temp_bear)

# interpolating lon and lat

from scipy.interpolate import interp1d

# goes through each key / range
for n, t, b, r in zip(temp_lon, temp_lat, temp_bear, range(rows)):

    # start of cluster
    start = 0
    
    # going through each bearing
    for i in range(1, len(temp_bear[b])):
    
        # checking if next point is too far or at the end
        if temp_bear[b][i] - temp_bear[b][i - 1] > 12 or i == len(temp_bear[b]) - 1:
            if i == len(temp_bear[b]) - 1:
                i = i + 1
            cluster = temp_bear[b][start:i]
            
            # checking if cluster has enough points
            if len(cluster) > 2:
                # +2 to finish the interpolation at the end
                new_angles = np.arange(cluster[0], cluster[-1] + 2, 2)
                if_lon = interp1d(np.array(cluster), np.array(temp_lon[n][start:i]), kind='quadratic')
                if_lat = interp1d(np.array(cluster), np.array(temp_lat[t][start:i]), kind='quadratic')
                i_lon = if_lon(new_angles)
                i_lat = if_lat(new_angles)

                # adding each interpolated lon and lat to the matrix
                for o, a, p in zip(i_lon, i_lat, new_angles):
                    
                    # getting the index for bearing
                    c = int((p - minc) / 2)
            
                    if interpolated_matrix[r][c] != 1:
                        lon_matrix[r][c] = o
                        lat_matrix[r][c] = a
                        interpolated_matrix[r][c] = 2
                
            # updating start of cluster
            start = i
         
# check for row or column interpolation (have neighbors)   
def check_next(r, c, l1, l2):
    if (interpolated_matrix[r-1][c] == 1 and interpolated_matrix[r+1][c] == 1 and interpolated_matrix[r][c-1] == 1 and interpolated_matrix[r][c+1]==1):
        l1.append(u_matrix[r-1][c])
        l1.append(u_matrix[r+1][c])
        l1.append(u_matrix[r][c-1])
        l1.append(u_matrix[r][c+1])
        l2.append(v_matrix[r-1][c])
        l2.append(v_matrix[r+1][c])
        l2.append(v_matrix[r][c-1])
        l2.append(v_matrix[r][c+1])
        return True
    elif (interpolated_matrix[r][c-1] == 1 and interpolated_matrix[r][c+1]==1):
        l1.append(u_matrix[r][c-1])
        l1.append(u_matrix[r][c+1])
        l2.append(v_matrix[r][c-1])
        l2.append(v_matrix[r][c+1])
        return True
    elif (interpolated_matrix[r-1][c] == 1 and interpolated_matrix[r+1][c] == 1):
        l1.append(u_matrix[r-1][c])
        l1.append(u_matrix[r+1][c])
        l2.append(v_matrix[r-1][c])
        l2.append(v_matrix[r+1][c])
        return True
    return False


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

# checking if grid is there
def check_grid(r, c, a):
    
    if (interpolated_matrix[r + a[0]][c - a[2]] != 1 or 
        interpolated_matrix[r - a[1]][c - a[2]] != 1 or
        interpolated_matrix[r + a[0]][c + a[3]] != 1 or
        interpolated_matrix[r - a[1]][c + a[3]] != 1):
            return False
    
    return True

# possible grids
def get_grid(r, c, a):
    
    # check [1, 1, 1, 1]
    if check_grid(r, c, a):
        return True

    # check when one entry is 2
    for e in a:
        e = 2
        if check_grid(r, c, a):
            return True
        e = 1
    
    # when two entries are 2
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

# bilinear interpolation formula
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
      
def add_values(r, c):
    interpolated_matrix[r][c] = 3
    lon_bilin.append(lon_matrix[r][c])
    lat_bilin.append(lat_matrix[r][c])
    
def add_uv(r, c, u, v):
    u_matrix[r][c] = u
    v_matrix[r][c] = v
    u_lin.append(u)
    v_lin.append(v)
    interpolated_matrix[r][c] = 3
    lon_lin.append(lon_matrix[r][c])
    lat_lin.append(lat_matrix[r][c])

# for average
from statistics import mean

# applying the interpolation
def interpolation(d):

    for r in range(d, rows - d):
    
        for c in range(d, cols - d):
            
            a = [1, 1, 1, 1]
            
            # if lon + lat has been interpolated
            if interpolated_matrix[r][c] == 2:
                
                l1 = []
                l2 = []
            
                # for points that are close enough
                if (check_next(r, c, l1, l2)):
                    u_result = mean(l1)
                    v_result = mean(l2)
                    add_uv(r, c, u_result, v_result)
            
                # add points if grid has been found
                elif get_grid(r, c, a):
                    bilinear(u_matrix, u_bilin, r, c, a)
                    bilinear(v_matrix, v_bilin, r, c, a)
                    add_values(r, c)

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
    #figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

# import plotting stuff
from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm, Normalize

# quiver configurations
scale=None
headwidth=3
headlength=5
headaxislength=5
velocity_min = -40
velocity_max = 40

# plot details
title = 'Radial Map: Original and Interpolated Data'
plt.title(f'{title}\n')

offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# turning into numpy arrays (combine both)
x = np.concatenate((lon_original, lon_lin, lon_bilin)) 
y = np.concatenate((lat_original, lat_lin, lat_bilin)) 
u = np.concatenate((u_original, u_lin, u_bilin))
v = np.concatenate((v_original, v_lin, v_bilin)) 

# setting the extent
#extent = [-82.8, -78.3, 22.5, 26.1]
extent = [x.min()-0.2, x.max()+0.2, y.min()-0.2, y.max()+0.2]
ax.set_extent(extent)

qargs = dict(scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# add color
colors = np.concatenate([['c'] * len(u_original), ['tomato'] * len(u_lin), ['b'] * len(u_bilin)])

# add legend
import matplotlib.lines as mlines
legend_lines = [mlines.Line2D([], [], color='c', lw=3, label='Original Data'), 
                 mlines.Line2D([], [], color='tomato', lw=3, label='Interpolated Data'),
                 mlines.Line2D([], [], color='b', lw=3, label='Bilinear Data')]

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

