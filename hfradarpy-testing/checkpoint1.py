from hfradarpy.radials import Radial, qc_radial_file
import glob
import os
import xarray as xr

# Path to radial directory
#radial_dir = '/home/cqiao/HFR_proc/hfradarpy/examples/data/radials/ruv/SEAB/'
#save_dir = '/home/cqiao/HFR_proc/hfradarpy/examples/data/radials_qc/ruv/SEAB'

#radial_dir = '/home/cqiao/HFR_proc/radials/'
radial_dir = './radials_raw/WEST/'
#save_dir = './radials_qc/nc/'

# Use glob to find radial files (*
#files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))
# files[:10]


radial_file = '/home/cqiao/HFR_data/radial_files/NODC_Test/RDL_i_USF_RDSR_2011_01_22_1800.hfrss10lluv'

#r = Radial(files[25])
r = Radial(radial_file)
r

# run high frequency radar qartod tests on open radial file

qc_values = dict(
    qc_qartod_avg_radial_bearing=dict(reference_bearing=151, warning_threshold=15, failure_threshold=30),
    qc_qartod_radial_count=dict(min_count=75.0, low_count=225.0),
    qc_qartod_maximum_velocity=dict(max_speed=300.0, high_speed=100.0),
    qc_qartod_spatial_median=dict(smed_range_cell_limit=2.1, smed_angular_limit=10, smed_current_difference=30),
    qc_qartod_temporal_gradient=dict(gradient_temp_fail=32, gradient_temp_warn=25),
    qc_qartod_primary_flag=dict(include=['qc_qartod_syntax', 'qc_qartod_valid_location', 'qc_qartod_radial_count',
                                         'qc_qartod_maximum_velocity', 'qc_qartod_spatial_median'])
)
r.initialize_qc()
r.qc_qartod_syntax()
r.qc_qartod_maximum_velocity(**qc_values['qc_qartod_maximum_velocity'])
r.qc_qartod_valid_location()
r.qc_qartod_radial_count(**qc_values['qc_qartod_radial_count'])
r.qc_qartod_spatial_median(**qc_values['qc_qartod_spatial_median'])
#r.qc_qartod_temporal_gradient(files[1]) #pass the previous hourly radial to this one
r.qc_qartod_avg_radial_bearing(**qc_values['qc_qartod_avg_radial_bearing'])
r.qc_qartod_primary_flag(**qc_values['qc_qartod_primary_flag'])

tds = r.to_xarray('gridded', enhance=True).squeeze()
tds

'''
save_file = save_dir + r.file_name[:-4:] + '.nc'
print('Saving NetCDF file: ' + save_file)
tds.to_netcdf(save_file)
'''

# Lets get rid of the single time dimension. It will cause problems during plotting
tds = tds.squeeze()

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

extent = []

dx = dy = 0.2  # Area around the point of interest.

extent = [tds.lon.min()-dx, tds.lon.max()+dx, tds.lat.min()-dy, tds.lat.max()+dy]

print(extent)

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

plt.title(f'Radial Plot\n{tds.time.data}')
plt.quiver(tds.lon.data, tds.lat.data, tds.u.data, tds.v.data, transform=ccrs.PlateCarree())

# Get the receiver location for plotting purposes
receiver_location = [float(x) for x in tds.Origin.split('  ')]
receiver_location.reverse()
receiver_location

plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())
#ax.set_extent([-83.0259, -78.9407, 22.8547, 26.6224])
map_features(ax)
plt.show()

import numpy.ma as ma

time = tds.time
lon = tds.coords['lon'].data
lat = tds.coords['lat'].data
u = tds['u'].data
v = tds['v'].data

u = ma.masked_invalid(u)
v = ma.masked_invalid(v)

import numpy as np
from oceans.ocfis import uv2spdir, spdir2uv

angle, speed = uv2spdir(u, v)  # convert u/v to angle and speed

u, v = spdir2uv(  # convert angle and speed back to u/v, normalizing the arrow sizes
    np.ones_like(speed),
    angle,
    deg=True
    )

# Intialize an empty subplot using cartopy
fig, ax = plt.subplots(
    figsize=(11, 8),
    subplot_kw=dict(projection=ccrs.Mercator())
)

plt.title(f'Normalized arrow sizes\n{tds.time.data}')
plt.quiver(lon, lat, u, v, transform=ccrs.PlateCarree())
plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())

map_features(ax)
plt.show()

import cmocean
from matplotlib.colors import TwoSlopeNorm, Normalize

"""
Displays the direction and magnitude of the radials
"""
cmap = cmocean.cm.balance
scale=50
headwidth=2.5
headlength=4
headaxislength=4
sub=1
velocity_min = -40
velocity_max = 40
cbar_step = 10
offset = Normalize(vmin=velocity_min, vmax=velocity_max, clip=True)

# Define arrow colors. Limited by velocity_min and velocity_max
color_clipped = np.clip(
    tds.velocity.data[::sub],
    velocity_min,
    velocity_max
).squeeze()

ticks = np.append(np.arange(velocity_min, velocity_max, cbar_step), velocity_max)

from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )

# Plot title
plt.title(f'Colorized Radial Plot\n{tds.time.data}')

qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# plot arrows over pcolor
h = ax.quiver(
    lon[::sub],
    lat[::sub],
    u[::sub],
    v[::sub],
    color_clipped,
    **qargs
)
map_features(ax)

# generate colorbar
divider = make_axes_locatable(ax)
cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
fig.add_axes(cax)

cb = plt.colorbar(h, cax=cax, ticks=ticks)
cb.ax.set_yticklabels([f'{s:d}' for s in ticks])
cb.set_label('cm/s')

plt.show()

from matplotlib import colors

"""
Motion displays the direction (towards or away) from radar
"""
title = 'Radial Map: Towards/Away from radar'
cmap= 'bwr'

velocity = tds.velocity
velocity_temp = velocity.where(velocity > 0, other=-1)  # Going away from radar
color_clipped = velocity_temp.where(velocity < 0, other=1).data  # Going towards radar
offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )

# Plot title
plt.title(f'{title}\n{tds.time.data}')

qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# plot arrows over pcolor
h = ax.quiver(
    lon[::sub],
    lat[::sub],
    u[::sub],
    v[::sub],
    color_clipped,
    **qargs
)
map_features(ax)
plt.show()

color_clipped = tds.primary_flag_qc.where(tds.primary_flag_qc == 1, other=-1).data  # PRIM == 1 where vectors pass qc
offset = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
title = "Radial Map: QC Pass/Fail"
cmap = colors.ListedColormap(['red', 'limegreen'])

fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )

# Plot title
plt.title(f'{title}\n{tds.time.data}')


qargs = dict(cmap=cmap, scale=scale, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength)
qargs['transform'] = ccrs.PlateCarree()
qargs['norm'] = offset

# plot arrows over pcolor
h = ax.quiver(
    lon[::sub],
    lat[::sub],
    u[::sub],
    v[::sub],
    color_clipped,
    **qargs
)
map_features(ax)
plt.colorbar(h)
plt.show()

