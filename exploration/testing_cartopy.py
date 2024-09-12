
# Created Jul 2024 by Cathleen Qiao

# This file is for testing cartopy functions as I kept getting a segmentation
# fault when trying to use ax.coastlines or ax.set_extent. I later installed
# miniconda onto my machine and everything worked fine. I then used this file
# to get a better understanding of what the functions could provide. 



import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs

#fig = plt.figure(figsize=(10, 5))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# make the map global rather than have it zoom in to
# the extents of any plotted data
# ax.set_global()

central_lon, central_lat = 100, -8.5

extent = [-50, 50, -30, 30]

#proj=ccrs.PlateCarree(central_longitude=100)

ax = plt.axes(projection=ccrs.PlateCarree())

#ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=100))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_lon, central_lat))

ax.set_extent(extent)
#ax.gridlines()

#ax.stock_img()
# add some features
ax.coastlines(resolution='50m')
#ax.coastlines(resolution='110m', color='black', linewidth=0.8)
ax.add_feature(cfeature.LAND)
#ax.add_feature(cfeature.OCEAN)
#ax.add_feature(cfeature.COASTLINE)
#ax.add_feature(cfeature.LAKES, alpha=0.5)
#ax.add_feature(cfeature.RIVERS)
#ax.add_feature(cfeature.BORDERS)

#ax.set_xlim(0,10)
#ax.add_feature(cfeature.LAND, edgecolor='black', facecolor="none", linewidth=1)

plt.show()
#plt.savefig('./test.jpg')
