from hfradarpy.radials import Radial, qc_radial_file
import glob
import os
import xarray as xr

# Path to radial directory
radial_dir = '/home/cqiao/HFR_proc/radials_raw/'
save_dir = '/home/cqiao/HFR_proc/pre_proc_images/radial_plots/'

#radial_dir = '/home/cqiao/HFR_proc/radials_qc/'
#save_dir = '/home/cqiao/HFR_proc/post_proc_images/radial_plots/'

# Use glob to find radial files (*
files = sorted(glob.glob(os.path.join(radial_dir, '*.ruv')))


def save_radial(r):
    tds = r.to_xarray('gridded', enhance=True).squeeze()
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

    extent = [tds.lon.min()-dx, tds.lon.max()+dx, tds.lat.min()-dy, tds.lat.max()-1]

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

    plt.title(f'Raw Data Radial Plot\n{tds.time.data}')
    #plt.title(f'Processed Data Radial Plot\n{tds.time.data}')
    plt.quiver(tds.lon.data, tds.lat.data, tds.u.data, tds.v.data, transform=ccrs.PlateCarree())

    # Get the receiver location for plotting purposes
    receiver_location = [float(x) for x in tds.Origin.split('  ')]
    receiver_location.reverse()
    receiver_location

    plt.plot(receiver_location[0], receiver_location[1], 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())
    map_features(ax)

    file_name = r.file_name[:-20:] + '_raw.png'
    #file_name = r.file_name[:-4:] + '.png'

    print(f'{file_name}')

    fig.savefig(save_dir + file_name)
    
    
    
for f in files[:1]:
    r = Radial(f)
    save_radial(r)

