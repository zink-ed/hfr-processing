# import packages
import numpy as np
import math
from pyproj import Geod
from shapely.geometry import Point
from geopandas import GeoSeries

'''
FUNCTIONS IN THIS FILE
    createGrid(lonMin, lonMax, latMin, latMax, gridRes)
    convertAngle(trueAngle, radians=False)
    dms2dd(dms)
    evaluateGDOP(cell, siteLon, siteLat, g)
'''

# createLonLatGridFromBB - creating grid of lon / lat coordinates ##############
def createGrid(lonMin, lonMax, latMin, latMax, gridRes):
    
    """ 
    this function creates a regular lon/lat grid given the limits
    for the bounding box and grid resolution (uses WGS84 CRS)
    
    INPUT: lonMin (westernmost)
           lonMax (easternmost)
           latMin (southermost)
           latMax (northernmost)
           gridRes (in meters)
    
    OUTPUT: pts = Geopandas GeoSeries containing lon / lat pairs of
                  all the points in the grid
            length = longitude range of grid
            width = latitude range of grid
    """
    
    # managing antimeridian crossing
    antiMeridianCrossing = False
    if lonMin > 0 and lonMax < 0:
        antiMeridianCrossing = True
        lonMax = lonMax + 360
        
    
    # use WGS84 ellipsoid
    g = Geod(ellps='WGS84')
    
    # forward and back azimuths, distance between lower left and upper left
    az12, az21, dist = g.inv(lonMin, latMin, lonMin, latMax)
    
    # retrieve array of distances from lower left points along lat axis
    dd = np.arange(0, dist, gridRes)
    
    # compute lat, lon, back azimuth of all points along lat axis
    fooLon, Lat, backaz = g.fwd(len(dd) * [lonMin], len(dd) * [latMin], len(dd) * [0], dd)
    Lat = np.array(Lat)
    
    # retrieve coordinates of center
    lonCenter = (lonMin + lonMax) / 2
    latCenter = (latMin + latMax) / 2
    
    # evaluate distance
    az12, az21, distLon = g.inv(lonCenter, latCenter, lonCenter+1e-4, latCenter)
    distLon = 1e4 * distLon
    
    # evaluate displacement in longitude
    dd = gridRes / distLon
    
    # computer lon of all points along lon axis
    Lon = np.arange(lonMin, lonMax, dd)
    
    # adjust lon and lat so they are centered in box
    Lon = Lon + (lonMax - Lon[-1]) / 2;
    Lat = Lat + (latMax - Lat[-1]) / 2;
    
    # manage antimeridian crossing
    if antiMeridianCrossing:
        Lon[Lon > 180] = Lon[Lon > 180] - 360
        
    # create grid
    length = len(Lon)
    width = len(Lat)
    #print(length)
    #print(width)
    Lon, Lat = np.meshgrid(Lon, Lat);
    Lonc = Lon.flatten()
    Latc = Lat.flatten()
    
    # convert these points to geo-data
    pts = GeoSeries([Point(x, y) for x, y in zip(Lonc, Latc)])
    pts = pts.set_crs('epsg:4326')
    
    return pts, length, width
    

# true2mathAngle - convert angles from true to math ############################
def convertAngle(trueAngle, radians=False):
    
    """
    this function converts angles from the true system (the geographicalsite_source
    system) to the math system (trigonometry system)
    
    math convention - an angle is measured CCW from East
    true convention - an angle is measured CW from North
    
    INPUT: trueAngle = numpy array containing true angles
           radians = flag indicating whether the angles are in degrees (default) 
                     or in radians
        
    OUTPUT: mathAngle = numpy array containing math angles
    """
    
    # convert to degrees if angles are in radians
    if radians:
        trueAngle = np.rad2deg(trueAngle)
        
    mathAngle = 90 - trueAngle
    mathAngle = np.mod(mathAngle,360)
    
    return mathAngle


def dms2dd(dms):
    
    """
    this function converts angles from DMS (degrees-minutes-seconds) to DD 
    (deciml degrees)
    
    INPUT: dms = tuple containing degree, minute and second values
        
    OUTPUT: dd = decimal degrees angle
    """
    
    dd = float(dms[0]) + float(dms[1])/60 + float(dms[2])/(60*60)
    
    return dd


def evaluateGDOP(cell, siteLon, siteLat, g):
    
    """
    this function evaluates the GDOP value of a grid cell based on its coordinates 
    and on the coordinates of the radial sites
    
    INPUT: cell = Series containing longitude and latitude of the grid cell for which the GDOP is evaluated
           siteLon = list containing the longitudes of the radial sites
           siteLat = list containing the latitudes of the radial sites
           g = Geod object with CRS.
        
    OUTPUT: gdop = GDOP value
    """
    
    # convert grid cell Series to numpy arrays
    cell = cell.to_numpy()
    cellLon = cell[0]
    cellLat = cell[1]
    
    # evaluate the radial angles from the radial sites to the grid cell
    radialAngles,az21,dist = g.inv(siteLon,siteLat,len(siteLon)*[cellLon],len(siteLat)*[cellLat])
    
    # form the design matrix for GDOP evaluation
    Agdop = np.stack((np.array([np.cos(np.deg2rad(radialAngles))]),np.array([np.sin(np.deg2rad(radialAngles))])),axis=-1)[0,:,:]
    
    # evaluate the covariance matrix Cgdop for GDOP evaluation
    Cgdop = np.linalg.inv(np.matmul(Agdop.T, Agdop))
    
    # evaluate GDOP
    gdop = math.sqrt(Cgdop.trace())
    
    return gdop



