# import libraries
from hfradarpy.radials import Radial
import glob
import os
import xarray as xr
import numpy as np

# read in the data of a radial file
radial_dir = './radials_clean/MARA/'

# use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

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

temp_lon = {}
temp_lat = {}
temp_bear = {}
temp_u = {}
temp_v = {}

def create_dict(list1, dict1):
    for r, l in zip(ranges, list1):
        if r in dict1:
            dict1[r].append(l)
        else:
            dict1[r] = [l]
        
create_dict(lon_original, temp_lon)
create_dict(lat_original, temp_lat)
#create_dict(lon_original, temp_lon)
#create_dict(lon_original, temp_lon)

i = 0

while ranges[i] == ranges[0]:
    #temp_lon.append(lon_original[i])
    temp_lat.append(lat_original[i])
    temp_bear.append(bearings[i])
    temp_u.append(u_original[i])
    temp_v.append(v_original[i])
    i = i + 1
   
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d, interp2d

import matplotlib.pyplot as plt


fig, ax = plt.subplots(
    figsize=(11, 8)
)


new_angles = np.arange(temp_bear[0], temp_bear[-1], 2)

# cubic spline and interp1d
'''
cf_lon = CubicSpline(np.array(temp_bear), np.array(temp_lon))
if_lon = interp1d(np.array(temp_bear), np.array(temp_lon), kind='quadratic')
c_lon = cf_lon(new_angles)
i_lon = if_lon(new_angles)

cf_lat = CubicSpline(np.array(temp_bear), np.array(temp_lat))
if_lat = interp1d(np.array(temp_bear), np.array(temp_lat), kind='quadratic')
c_lat = cf_lat(new_angles)
i_lat = if_lat(new_angles)
'''
# interp2d

'''
if_u = interp2d(np.array(temp_bear), np.array(temp_u), kind='quadratic')
if_v = interp2d(np.array(temp_bear), np.array(temp_v), kind='quadratic')
i_u = if_u(new_angles)
i_v = if_v(new_angles)
'''

'''
# try rbf interpolation
r = np.deg2rad(temp_bear)


# plotting
plt.plot(new_angles, i_u, color = 'r', alpha=0.75)
plt.plot(new_angles, c_u, color ='g', alpha=0.65)
plt.scatter(temp_bear, temp_u, color ='b')
plt.show()
'''

'''
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
    
'''
