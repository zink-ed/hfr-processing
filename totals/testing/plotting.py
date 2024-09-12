
# Created Aug 2024 by Cathleen Qiao

# This file is for plotting totals using cartopy. I wanted to implement this to
# see which plotting method would look nicer. Additionally, if I was able to 
# plot the totals with cartopy, they would match the graphs of the radials. 

# This file is unfinished. 


# import plotting packages
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    
         
def plot_cartopy(T):

    # set colors of the land. 
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
    
    # Intialize an empty subplot using cartopy
    fig, ax = plt.subplots(
        figsize=(11, 8),
        subplot_kw=dict(projection=ccrs.Mercator())
    )
    
    # Compute the native map projection coordinates for the vectors
    x, y = (T.data.LOND, T.data.LATD)
    
    extent = [x.min()-0.2, x.max()+0.2, y.min()-0.2, y.max()+0.2]
    ax.set_extent(extent)
    
    # Create the velocity component variables
    u = T.data.VELU / 100        # CODAR velocities are in cm/s
    v = T.data.VELV / 100        # CODAR velocities are in cm/s 
    vel = abs(T.data.VELO) / 100  

    plt.title(f'Totals')
    plt.quiver(x, y, u, v, vel, cmap=plt.cm.jet, width=0.001, headwidth=4, headlength=4, headaxislength=4)

    # get station coordinates and codes
    siteLon = T.site_source['Lon'].to_numpy()
    siteLat = T.site_source['Lat'].to_numpy()
    siteCode = T.site_source['Name'].tolist()
    
    # compute the native map projection coordinates for the stations
    xS, yS = (siteLon, siteLat)       
    
    # plot radial stations
    for label, xs, ys in zip(siteCode, xS, yS):
        plt.plot(xs, ys, 'o', markersize=10, markeredgecolor='black', color='red', transform=ccrs.PlateCarree())
        ax.annotate(label, (xs, ys))


    # axes properties and features
    ax.add_feature(LAND, edgecolor=edgecolor, facecolor=landcolor)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, zorder=11, edgecolor=edgecolor)

    # gridlines and grid labels
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
    
    plt.show()

