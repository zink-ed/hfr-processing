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
create_dict(bearings, temp_bear)
#create_dict(lon_original, temp_lon)

#print(temp_lon)
#print(temp_lat)
#print(temp_bear)

'''
i = 0

while ranges[i] == ranges[0]:
    #temp_lon.append(lon_original[i])
    temp_lat.append(lat_original[i])
    temp_bear.append(bearings[i])
    temp_u.append(u_original[i])
    temp_v.append(v_original[i])
    i = i + 1
'''

from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d, interp2d

import matplotlib.pyplot as plt


fig, ax = plt.subplots(
    figsize=(11, 8)
)


#new_angles = np.arange(temp_bear[0], temp_bear[-1], 2)

# cubic spline and interp1d


# goes through each key
for n, t, b in zip(temp_lon, temp_lat, temp_bear):
    i_lon = []
    i_lat = []
    
    if len(temp_bear[b]) > 25:
        #print(len(temp_bear[b]))
        new_angles = np.arange(temp_bear[b][0], temp_bear[b][-1], 2)
        if_lon = interp1d(np.array(temp_bear[b]), np.array(temp_lon[n]), kind='quadratic')
        if_lat = interp1d(np.array(temp_bear[b]), np.array(temp_lat[t]), kind='quadratic')
        i_lon = if_lon(new_angles)
        i_lat = if_lat(new_angles)
        plt.scatter(i_lon, i_lat, color='r', alpha=0.75)
    plt.scatter(temp_lon[n], temp_lat[t], color='b', alpha=0.25)
    

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
'''

# plotting
#plt.plot(new_angles, i_u, color = 'r', alpha=0.75)
#plt.plot(new_angles, c_u, color ='g', alpha=0.65)
#plt.scatter(temp_bear, temp_u, color ='b')
plt.show()


