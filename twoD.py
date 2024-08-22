# import libraries
from hfradarpy.radials import Radial
import glob
import os
import xarray as xr

# read in the data of a radial file
radial_dir = './radials_clean/MARA/'
#radial_dir = './radials_clean/KEYWEST/'

# use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))

# select file
r = Radial(files[70])

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

# interpolating lon and lat

from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d, interp2d, griddata

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
            
            #print(cluster)
            
            # checking if cluster has enough points
            if len(cluster) > 2:
                # +2 to finish the interpolation at the end
                new_angles = np.arange(cluster[0], cluster[-1] + 2, 2)
                if_lon = interp1d(np.array(cluster), np.array(temp_lon[n][start:i]), kind='quadratic')
                if_lat = interp1d(np.array(cluster), np.array(temp_lat[t][start:i]), kind='quadratic')
                i_lon = if_lon(new_angles)
                i_lat = if_lat(new_angles)
                #plt.scatter(i_lon, i_lat, color='r', alpha=0.75)

                # adding each interpolated lon and lat to the matrix
                for o, a, p in zip(i_lon, i_lat, new_angles):
                    
                    # getting the index for bearing
                    c = int((p - minc) / 2)
            
                    if interpolated_matrix[r][c] != 1:
                        lon_matrix[r][c] = o
                        lat_matrix[r][c] = a
                        lon_interpolated.append(o)
                        lat_interpolated.append(a)
                        interpolated_matrix[r][c] = 2
                
            # updating start of cluster
            start = i
                        
    #plt.scatter(temp_lon[n], temp_lat[t], color='b', alpha=0.25)

x = np.concatenate((lon_original, lon_interpolated)) 
y = np.concatenate((lat_original, lat_interpolated)) 

#print(len(lon_original))
#print(len(lon_interpolated))

points=np.transpose(np.vstack((lon_original, lat_original)))
grid = np.transpose(np.vstack((x, y)))

#print(points)
#print(grid)

u_m2 = griddata(points, u_original, grid, method='cubic')
v_m2 = griddata(points, v_original, grid, method='cubic')
#print(u_original)
#print(v_original)
#print(u_m2)
#print(v_m2)

# plotting from the matrix (to monitor progress)
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
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

offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

# turning into numpy arrays (combine both)
#u = np.concatenate((u_original, u_interpolated))
#v = np.concatenate((v_original, v_interpolated)) 

# setting the extent
extent = [-82.8, -78.3, 22.5, 26.1]
#extent = [x.min()-0.2, x.max()+0.2, y.min()-0.2, y.max()+0.2]
ax.set_extent(extent)

qargs = dict(scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# add color
colors = np.concatenate([['c'] * len(lon_original), ['tomato'] * (len(lon_interpolated))])

# plot arrows
q1 = ax.quiver(
    x,
    y,
    u_m2,
    v_m2,
    color=colors,
    **qargs
)

#plt.plot(new_angles, i_u, color = 'r', alpha=0.75)
#plt.plot(new_angles, c_u, color ='g', alpha=0.65)
#plt.scatter(x, y, color ='r', alpha=0.65)
#plt.scatter(lon_original, lat_original, color ='b', alpha=0.45)
plt.show()
