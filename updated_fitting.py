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
r = Radial(files[10])

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
    interpolated_matrix[r][b] = 1
    
    
# getting interpolated data

u_interpolated = []
v_interpolated = []

temp_lon = {}
temp_lat = {}
temp_bear = {}


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

# interpolated lon and lat

lon_interpolated = []
lat_interpolated = []

from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d, interp2d

# plotting from the matrix (to monitor progress)
import matplotlib.pyplot as plt


# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    #subplot_kw=dict(projection=ccrs.Mercator())
)

# goes through each key
for n, t, b in zip(temp_lon, temp_lat, temp_bear):

    
    i_lon = []
    i_lat = []
    
    start = 0
    print("NEW")
    print(temp_bear[b])
    
    for i in range(1, len(temp_bear[b])):
        if temp_bear[b][i] - temp_bear[b][i - 1] > 12 or i == len(temp_bear[b]) - 1:
            if i == len(temp_bear[b]) - 1:
                i = i + 1
            cluster = temp_bear[b][start:i]
            
            print(cluster)
            
            if len(cluster) > 2:
                #print(len(temp_bear[b]))
                new_angles = np.arange(cluster[0], cluster[-1] + 2, 2)
                if_lon = interp1d(np.array(cluster), np.array(temp_lon[n][start:i]), kind='quadratic')
                if_lat = interp1d(np.array(cluster), np.array(temp_lat[t][start:i]), kind='quadratic')
                i_lon = if_lon(new_angles)
                i_lat = if_lat(new_angles)
                plt.scatter(i_lon, i_lat, color='r', alpha=0.75)
                
            start = i
            
    plt.scatter(temp_lon[n], temp_lat[t], color='b', alpha=0.25)


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

#split()

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
#interpolation(2)    

# turning into numpy arrays (combine both)
x = np.concatenate((lon_original, lon_interpolated)) 
y = np.concatenate((lat_original, lat_interpolated)) 
u = np.concatenate((u_original, u_interpolated))
v = np.concatenate((v_original, v_interpolated)) 


# include legend
#plt.scatter(x, y)
#ax.legend(handles=legend_lines)
plt.show()

