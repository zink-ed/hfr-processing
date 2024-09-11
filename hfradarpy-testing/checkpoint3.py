from hfradarpy.radials import Radial, concat
import glob
import os
import xarray as xr

# Path to radial directory
radial_dir = '/home/cqiao/HFR_proc/radials_raw/'

# Use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))
files[:10] # List first 10 radials found for brevity

ds = concat(sorted(files), method='gridded', enhance=True)
ds
# ds = xr.open_mfdataset(radial_dir + '*.nc')

# Let's plot the latest hour available
tds = ds.isel(time=-1) # use isel instead of sel to select the time by index rather than the time you want. In this case, -1 will grab the last time
tds

# Get the receiver location for plotting purposes
receiver_location = [float(x) for x in ds.Origin.split('  ')]
receiver_location.reverse()
receiver_location

# Import matplotlib and cartopy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Adjust some standard plotting settings to make them the size of a sheet of paper
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 8
plt.rcParams["figure.figsize"] = fig_size

extent = []

dx = dy = 0.2  # Area around the point of interest.

extent = [ds.lon.min()-dx, ds.lon.max()+dx, ds.lat.min()-dy, ds.lat.max()+dy]

extent = [
    receiver_location[0] - dx, 
    receiver_location[0] + dx, 
    receiver_location[1] - dy, 
    receiver_location[1] + dy
]

print(extent)

# Set colors of the land. 
edgecolor = 'black'
landcolor = 'tan'

LAND = cfeature.NaturalEarthFeature(
    'physical', 'land', '10m',
    edgecolor='face',
    facecolor='tan'
)

state_lines = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
)

# Create a re-usable function for map features that we can pass an axes to.

def map_features(ax):
    # Axes properties and features
    ax.set_extent(extent)
    ax.add_feature(LAND, edgecolor=edgecolor, facecolor=landcolor)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, zorder=11, edgecolor=edgecolor)

    # Gridlines and grid labels
    gl = ax.gridlines(
        draw_labels=True,
        linewidth=.5,
        color='black',
        alpha=0.25,
        linestyle='--'
    )

    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}
    gl.xlocator = mticker.MaxNLocator(integer=True)
    gl.ylocator = mticker.MaxNLocator(integer=True)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.tick_params(which='major',
                   direction='out',
                   bottom=True, top=True,
                   labelbottom=True, labeltop=False,
                   left=True, right=True,
                   labelleft=True, labelright=False,
                   length=5, width=2)

    ax.tick_params(which='minor',
                   direction='out',
                   bottom=True, top=True,
                   labelbottom=True, labeltop=False,
                   left=True, right=True,
                   labelleft=True, labelright=False,
                   width=1)

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

plt.quiver(tds.lon.data, tds.lat.data, tds.u.data, tds.v.data, transform=ccrs.PlateCarree())
plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())

map_features(ax)

point_to_find = [-80.8, 24.6] #SEAB
#point_to_find = [-80.5, 24.5] #MARA

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

plt.quiver(tds.lon.data, tds.lat.data, tds.u.data, tds.v.data, transform=ccrs.PlateCarree())
plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())
plt.plot(point_to_find[0], point_to_find[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())

map_features(ax)

import pyproj

geodesic = pyproj.Geod(ellps='WGS84') #define the coordinate system. WGS84 is the standard used by GPS.

# Inverse transformation
# Determine forward and back azimuths, plus distances between initial points and terminus points.
_, back_azimuth, distance = geodesic.inv(receiver_location[0], receiver_location[1], point_to_find[0], point_to_find[1])

print(f'Bearing (CCWE): {back_azimuth} degrees, Range: {distance} meters') # degrees, meters

if back_azimuth < 0:
    back_azimuth = back_azimuth + 180
    
distance = distance/1000 #convert from meters to kilometers

print(f'Bearing (CWN): {back_azimuth} degrees, Range: {distance} kilometers') # degrees, meters

# Select nearest neighbor based on bearing and range
tds = ds.sel(bearing=back_azimuth, range=distance, method='nearest')
tds

# Plot of single grid point nearest to point selected above
tds.velocity.plot()
plt.grid()
# axes[0].set_ylim([-7, 7])
plt.xlabel('')
plt.ylabel('Velocity (cm/s)')
plt.suptitle(f'Radial Velocity at {point_to_find[1]}, {point_to_find[0]}')
plt.show()

# Use slice selection to get all points between selected bearings and ranges
bearing_slice = [115, 125]
range_slice = [50, 60]

tds = ds.sel(bearing=slice(bearing_slice[0], bearing_slice[1]), range=slice(range_slice[0], range_slice[1]))
tds

# Mean of points selected in previous cell
mean = tds.mean(dim=('bearing', 'range'))
mean

# Velocity Plot of mean
mean.velocity.plot()
plt.grid()
plt.xlabel('')
plt.ylabel('Velocity (cm/s)')
plt.suptitle(f'Average Radial Velocity\n All points between/including a bearing of {bearing_slice[0]} to {bearing_slice[1]} degrees and a range of {range_slice[0]} to {range_slice[1]} km ')
plt.show()

fig, axes = plt.subplots(nrows=2)
plt.suptitle(f'Average Radial Velocity\n All points between/including a bearing of {bearing_slice[0]} to {bearing_slice[1]} degrees and a range of {range_slice[0]} to {range_slice[1]} km ')

mean.u.plot(ax=axes[0])
axes[0].grid()
# axes[0].set_ylim([-7, 7])
axes[0].set_xlabel('')
axes[0].set_ylabel('u (cm/s)')

mean.v.plot(ax=axes[1])
axes[1].grid()
axes[1].set_ylim([-80, 80])
axes[1].set_ylabel('v (cm/s)')
axes[1].set_xlabel('')

plt.tight_layout()
plt.show()
# plt.savefig('/Users/mikesmith/Desktop/isaias-mean.png', dpi=300, bbox_inches='tight', pad_inches=0.1)

# Single bearing bin, multiple ranges
azimuth = 0
bearing = tds.isel(bearing=azimuth)

fig, axes = plt.subplots(nrows=2, sharex=True)
bearing.u.plot.line(x='time', ax=axes[0], )
axes[0].grid()
# axes[0].set_ylim([-7, 7])
axes[0].set_ylabel('u (cm/s)')
axes[0].set_xlabel('')

bearing.v.plot.line(x='time', ax=axes[1],)
axes[1].grid()
axes[1].set_ylim([-80, 80])
axes[1].set_ylabel('v (cm/s)')
axes[1].set_xlabel('')
plt.suptitle(f'Radial Velocity\n All points at a bearing of {azimuth} degrees and ranges including {bearing.range.data} km ')

plt.tight_layout()
plt.show()


