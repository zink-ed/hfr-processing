
from scipy.interpolate import griddata
import numpy as np
from calc import createGrid
from geopandas import GeoSeries
from shapely.geometry import Point    
from shapely.ops import nearest_points
   
 
def interpolation(T):
    
    lon_orig = T.data.LOND
    lat_orig = T.data.LATD
    
    points_orig = GeoSeries([Point(x, y) for x, y in zip(lon_orig, lat_orig)])
    #print(points_orig)
    
    points = np.transpose(np.vstack((lon_orig, lat_orig)))
    #print(points)
    
    
    # -----
    
    g, l, w = createGrid(min(lon_orig), max(lon_orig), min(lat_orig), max(lat_orig), 3000)
    #print(g)
    
    for p in g:
        #print(p)
        n = nearest_points(points_orig, p)
        print(n)
        print(n.shape)
        if n.x - p.x > 0.5 or n.y - p.y > 0.5:
            g.pop(p)
    
    x = np.unique(g.x.to_numpy())
    y = np.unique(g.y.to_numpy())
    
    lon, lat = np.meshgrid(x, y)
    lonc = lon.flatten()
    latc = lat.flatten()
    
    grid = np.transpose(np.vstack((lonc, latc)))
    #print(grid)
   

    u_orig = T.data.VELU / 100       # CODAR velocities are in cm/s
    v_orig = T.data.VELV / 100  
    vel_orig = T.data.VELO / 100
    
    u = griddata(points, u_orig, grid, method='cubic')
    v = griddata(points, v_orig, grid, method='cubic')
    vel = griddata(points, vel_orig, grid, method='cubic')
    
    return lonc, latc, u, v, vel

