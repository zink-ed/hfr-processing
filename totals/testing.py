import xarray as xr
import pandas as pd
import numpy as np

from radials import Radial
from totals import Total

# read in the data of a radial file
radial_dir1 = '../radials_clean/MARA/'
radial_dir2 = '../radials_clean/KEYWEST/'

import glob
import os

# search for radial files
files1 = sorted(glob.glob(os.path.join(radial_dir1, '*.ruv')))
files2 = sorted(glob.glob(os.path.join(radial_dir2, '*.ruv')))

r1 = Radial(radial_dir1 + 'RDLm_MARA_2024_04_16_1600.ruv')
r2 = Radial(radial_dir2 + 'RDLm_WEST_2024_04_16_1600.ruv')

#r1.to_xarray_multidimensional()
#r2.to_xarray_multidimensional()

#r1.plot()
#r2.plot()

# SELECT RADIALS

rDF = pd.DataFrame(columns=['Radial'])



# ----------------------------- FUNCTIONS --------------------------------------

from pyproj import Geod
from shapely.geometry import Point
from geopandas import GeoSeries

# createLonLatGridFromBB - creating grid of lon / lat coordinates --------------
def createGrid(lonMin, lonMax, latMin, latMax, gridRes):
    
    """ 
    this function creates a regular lon/lat grid given the limits
    for the bounding box and grid resolution (uses WGS84 CRS)
    
    INPUT: lonMin (westernmost)
           lonMax (easternmost)
           latMin (southermost)
           latMax (northernmost)
           gridRes (in meters)
    
    OUTPUT: pts - Geopandas GeoSeries containing lon / lat pairs of
                  all the points in the grid
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
    Lon, Lat = np.meshgrid(Lon, Lat);
    Lonc = Lon.flatten()
    Latc = Lat.flatten()
    
    # convert these points to geo-data
    pts = GeoSeries([Point(x, y) for x, y in zip(Lonc, Latc)])
    pts = pts.set_crs('epsg:4326')
    
    return pts
    
# testing function
lon_min = min(r1.data.LOND.min(), r2.data.LOND.min())
lon_max = max(r1.data.LOND.max(), r2.data.LOND.max())
lat_min = min(r1.data.LATD.min(), r2.data.LATD.min())
lat_max = max(r1.data.LATD.max(), r2.data.LATD.max())
gridRes = 3000
grid = createGrid(lon_min, lon_max, lat_min, lat_max, gridRes)
#print(grid)




# true2mathAngle - convert angles from true to math ----------------------------
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




'''
def buildTotal():

def createTotal():
    
    T = buildTotal()

'''


















################################################################################

import glob
import os

# selectRadials - create DataFrame of radial files to be combined --------------
def selectRadials(paths, time):
    
    """
    this function lists the input radial files that are within the processing
    time interval
    
    creates the DataFrame containing info needed for the combination of radial
    files into totals
    
    INPUTS: paths = list of possible paths of data
            time = datetime of the processing period (string)
 
    OUTPUTS: radialsTBP = DataFrame containing all the radials 
             to be processed for the input  
    """
    
    # create output total Series
    radialsTBP = pd.DataFrame(columns=['filename', 'filepath', 'timestamp'])
    
    files = []
    
    wildcard = '*' + time + '.ruv'
    
    for p in paths:
        # search for radial files
        files.extend((glob.glob(os.path.join(p, wildcard))))
    
    #print(files)
    
    # list all radial files                
    for f in files:
        
        # get file parts
        filePath = os.path.dirname(f)
        fileName = os.path.basename(f)
        
        # get file timestamp
        radial = Radial(f)
        timeStamp = radial.time.strftime("%Y %m %d %H %M %S")                    
        dateTime = radial.time.strftime("%Y-%m-%d %H:%M:%S")  

        # prepare data to be inserted into the output DataFrame
        dataRadial = {'filename': [fileName], 'filepath': [filePath], 'timestamp': [timeStamp], 'datetime': [dateTime]}
        
        dfRadial = pd.DataFrame(dataRadial)
        
        # insert into the output DataFrame
        radialsTBP = pd.concat([radialsTBP, dfRadial], ignore_index=True)
    
    return radialsTBP

# testing
rDF = selectRadials([radial_dir1, radial_dir2], '2024_04_16_1600')
print(rDF)

################################################################################

# processRadials - combine radials into totals
def processRadials(rads):
    
    """
    this function processes the input radial files and then combines them into totals
    
    INPUTS: rads = DataFrame containing the radials to be processed grouped by timestamp
                         for the input network with the related information

    OUTPUTS:
        
    """

    # add Radial objects to the DataFrame
    rads['Radial'] = (rads.filepath + '/' + rads.filename).apply(lambda x: Radial(x))

'''    
    # Rename indices with site codes
    indexMapper = dict(zip(rads.index.values.tolist(),rads['station_id'].to_list()))
    rads.rename(index=indexMapper,inplace=True)    
    
          
    # convert Radials to standard data format (netCDF)
    
        # European standard data model
        rads = rads.apply(lambda x: applyEHNradialDataModel(x, networkData,             stationData.loc[stationData['station_id'] == x.station_id],vers,logger), axis=1)
        
            
    # combine Radials into Total
    
    # check if at least two Radials are available for combination
    if len(groupedRad) > 1:
        dfTot = performRadialCombination(groupedRad,networkData,vers,logger)                
    
    if 'dfTot' in locals():
          
        # apply QC to Totals
        dfTot['Total'] = dfTot.apply(lambda x: applyEHNtotalQC(x, networkData, vers, logger), axis=1)        
            
    # convert Totals to standard data format (netCDF)
    
        # European standard data model
        dfTot = dfTot.apply(lambda x: applyEHNtotalDataModel(x, networkData, stationData, vers, logger), axis=1)

         
    return

'''

################################################################################

# RadBinsInSearchRadius
def findBins(cell, radial, sR, g):
    
    """
    this function finds out which radial bins are within the spatthresh of the
    origin grid cell (the WGS84 CRS is used for distance calculations)
    
    INPUT: cell = Series containing lons and lats of the origin grid cells
           radial = Radial object
           sR = search radius (in meters)
           g = Geod object according to the Total CRS
        
    OUTPUT: radInSr = list of the radial bins falling within the search radius of the
                      origin grid cell
    """
    
    # convert grid cell Series and radial bins DataFrame to numpy arrays
    cell = cell.to_numpy()
    radLon = radial.data['LOND'].to_numpy()
    radLat = radial.data['LATD'].to_numpy() 
    
    # evaluate distances between origin grid cells and radial bins
    az12,az21,cellToRadDist = g.inv(len(radLon) * [cell[0]],len(radLat) * [cell[1]], radLon, radLat)
    
    # figure out which radial bins are within the spatthresh of the origin grid cell
    radInSR = np.where(cellToRadDist < sR)[0].tolist()
    
    #if (len(radInSR) > 1):
    #    print(radInSR)
    
    return radInSR



# combineRadials
def combineRadials(rDF, grid, search, gRes, tStp):
    
    """
    this function generataes total vectors from radial measurements using the
    weighted Least Square method for combination (need at least 2 sites)
    
    INPUT: rDF = DataFrame containing input Radials; indices must be the site codes.
           grid = grid of all the lon / lat points
           search = search radius for combination in meters
           gRes = grid resoultion (in meters)
           tStp = timestamp in datetime format (YYYY-MM-DD hh:mm:ss)
        
    OUTPUT: Tcomb = Total object generated by the combination
    """

    # create empty total with grid
    Tcomb = Total(grid=grid)
    
    processRadials(rDF)

    # fill site_source DataFrame with contributing radials information
    siteNum = 0    # initialization of site number
    for index, row in rDF.iterrows():
        print(index)
        print(row)

        siteNum = siteNum + 1
        rad = row['Radial']
        thisRadial = pd.DataFrame(index=[index], columns=['#', 'Name', 'Lat', 'Lon', 'Coverage(s)', 'RngStep(km)', 'Pattern', 'AntBearing(NCW)'])
        thisRadial['#'] = siteNum
        thisRadial['Name'] = index
        thisRadial['Lat'] = float(rad.metadata['Origin'].split()[0])
        thisRadial['Lon'] = float(rad.metadata['Origin'].split()[1])
        thisRadial['Coverage(s)'] = float(rad.metadata['TimeCoverage'].split()[0])
        thisRadial['RngStep(km)'] = float(rad.metadata['RangeResolutionKMeters'].split()[0])
        thisRadial['Pattern'] = rad.metadata['PatternType'].split()[0]
        thisRadial['AntBearing(NCW)'] = float(rad.metadata['AntennaBearing'].split()[0])
    Tcomb.site_source = pd.concat([Tcomb.site_source, thisRadial])
            
    # insert timestamp
    Tcomb.time = tStp

    # fill Total with some metadata
    Tcomb.metadata['TimeZone'] = rad.metadata['TimeZone']
    Tcomb.metadata['AveragingRadius'] = str(search / 1000) + ' km'
    Tcomb.metadata['GridAxisOrientation'] = '0.0 DegNCW'
    Tcomb.metadata['GridSpacing'] = str(gRes / 1000) + ' km'

    # Create Geod object according to the Total CRS
    g = Geod(ellps=Tcomb.metadata['GreatCircle'].split()[0])

    # create DataFrame for storing indices of radial bins within the search radius of each cell
    combineRadBins = pd.DataFrame(columns=range(len(Tcomb.data.index)))


    # figure out which radial bins are within the spatthresh of each grid cell
    for index, row in rDF.iterrows():
        rad = row['Radial']         
        thisRadBins = Tcomb.data.loc[:,['LOND','LATD']].apply(lambda x: findBins(x, rad, search, g), axis=1)
        combineRadBins.loc[index] = thisRadBins
      
    #print(Tcomb)
    

    # loop over grid points and pull out contributing radial vectors
    combineRadBins = combineRadBins.T
    totData = combineRadBins.apply(lambda x: makeTotalVector(x,rDF), axis=1)

    # assign column names to the combination DataFrame
    totData.columns = ['VELU', 'VELV','VELO','HEAD','UQAL','VQAL','CQAL','GDOP','NRAD']

    # fill Total with combination results
    Tcomb.data[['VELU', 'VELV','VELO','HEAD','UQAL','VQAL','CQAL','GDOP','NRAD']] = totData
    
 
    # mask out vectors on land
    #Tcomb.mask_over_land(subset=True)
    
    # get the indexes of grid cells without total vectors
    indexNoVec = Tcomb.data[Tcomb.data['VELU'].isna()].index
    
    # delete these row indexes from DataFrame
    Tcomb.data.drop(indexNoVec, inplace=True)
    Tcomb.data.reset_index(level=None, drop=False, inplace=True)    
        # Set drop=True if the former indices are not necessary
        
    return Tcomb


# testing
#combineRadials(rDF, grid, gridRes)



# totalLeastSquare
def totalLeastSquare(velDF):
    
    """
    this function calculates the u/v components of a total vector from 2 to n 
    radial vector components using weighted Least Square method
    
    INPUT: velDF = DataFrame containing contributing radial velocities, bearings
                        and standard deviations
        
    OUTPUT: u = U component of the total vector
            v = V component of the total vector
            cov = covariance matrix
            covGDOP = covariance matrix assuming uniform unit errors 
                    for all radials (i.e. all radial std=1)
    """
    
    # convert angles from true convention to math convention
    velDF['HEAD'] = convertAngle(velDF['HEAD'].to_numpy())
    
    # form the design matrix (i.e. the angle matrix)
    A = np.stack((np.array([np.cos(np.deg2rad(velDF['HEAD'])) / velDF['STD']]), np.array([np.sin(np.deg2rad(velDF['HEAD'])) / velDF['STD']])), axis = -1)[0,:,:]
    
    # form the velocity vector
    b = (velDF['VELO'].to_numpy()) / velDF['STD']    
    
    # evaluate the covariance matrix cov (variance(U) = cov(1,1) and variance(V) = cov(2,2))
    A2 = np.matmul(A.T, A)
    if np.linalg.det(A2) > 0:
        cov = np.linalg.inv(A2)
    
        # calculate the u and v for the total vector
        a = np.matmul(cov, np.matmul(A.T, b))
        u = a[0]
        v = a[1]    
        
        # form the design matrix for GDOP evaluation (i.e. setting all radial std to 1)
        des = np.stack((np.array([np.cos(np.deg2rad(velDF['HEAD']))]), np.array([np.sin(np.deg2rad(velDF['HEAD']))])), axis = -1)[0,:,:]
        
        # evaluate the covariance matrix covGDOP for GDOP evaluation (i.e. setting all radial std to 1)
        des = np.matmul(des.T, des)
        if np.linalg.det(des):
            covGDOP = np.linalg.inv(des)
            #print(cov)
            #print(covGDOP)
            
            return u, v, cov, covGDOP
        
        else:
            u = np.nan
            v = np.nan
            cov = np.nan
            covGDOP = np.nan
            return u, v, C, covGDOP
    
    else:
        u = np.nan
        v = np.nan
        cov = np.nan
        covGDOP = np.nan
        return u, v, cov, covGDOP
        

import math

# makeTotalVector
def makeTotalVector(rBins, rDF):
    
    """
    this function gets the total vector for each grid cell
    (use weighted Least Square method for combination)
    
    INPUT: rBins = Series containing contributing radial indices
           rDF = DataFrame containing input Radials
        
    OUTPUT: totalData = Series containing u/v components and related errors of 
                        total vector for each grid cell
    """
    
    # set minimum number of contributing radial sites
    minContrSites = 2
    
    # set minimum number of contributing radial vectors
    minContrRads = 3
    
    # create output total Series
    totalData = pd.Series(np.nan,index=range(9))
    
    # only consider contributing radial sites
    contrRad = rBins[rBins.str.len() != 0]
    
    # check if there are at least two contributing radial sites
    if contrRad.size >= minContrSites:
    
        # loop over contributing radial indices for collecting velocities and angles
        contributions = pd.DataFrame()
        
        for idx in contrRad.index:
            contrVel = rDF.loc[idx]['Radial'].data.VELO[contrRad[idx]]
            contrHead = rDF.loc[idx]['Radial'].data.HEAD[contrRad[idx]]
            contrStd = rDF.loc[idx]['Radial'].data.ETMP[contrRad[idx]]
            contrStd = contrStd.rename("STD")
            contributions = pd.concat([contributions, pd.concat([contrVel, contrHead, contrStd], axis=1)])
                    
        # only keep contributing radials with valid standard deviation values
        contributions = contributions[contributions.STD.notnull()]
        contributions = contributions[contributions.STD != 0]
        
        # check if there are at least three contributing radial vectors
        if len(contributions.index) >= minContrRads:
            
            # combine radial contributions to get total vector for the current grid cell
            u, v, cov, covGDOP = totalLeastSquare(contributions)
            
            if not math.isnan(u):
                # populate Total Series
                totalData.loc[0] = u                                            # VELU
                totalData.loc[1] = v                                            # VELV
                totalData.loc[2] = np.sqrt(u**2 + v**2)                         # VELO
                totalData.loc[3] = (360 + np.arctan2(u,v) * 180/np.pi) % 360    # HEAD
                totalData.loc[4] = math.sqrt(cov[0,0])                          # UQAL
                totalData.loc[5] = math.sqrt(cov[1,1])                          # VQAL
                totalData.loc[6] = cov[0,1]                                     # CQAL
                totalData.loc[7] = math.sqrt(np.abs(covGDOP.trace()))           # GDOP
                totalData.loc[8] = len(contributions.index)                     # NRAD
            
    #print(totalData)
    return totalData












































import datetime as dt

bb = [lon_min, lon_max, lat_min, lat_max, gridRes]


# performRadialCombination - combine radials into totals
def performRadialCombination(combRad, bb):
    
    """
    this function performs the least square combination of the input Radials and creates
    a Total object containing the resulting total current data
    
    the function creates a DataFrame containing the resulting Total object along with 
    related information
    
    INPUTS: combRad = DataFrame containing the Radial objects to be combined
            bb = list containing info about bounding box and grid resolution

    OUTPUTS: combTot = DataFrame containing the Total object obtained
        
    """
    
    # create the output DataFrame
    combTot = pd.DataFrame(columns=['Total'])
    
    # create the geographical grid
    grid = createGrid(bb[0], bb[1], bb[2], bb[3], bb[4])
    
    # get the combination search radius in meters
    #searchRadius = networkData.iloc[0]['combination_search_radius'] * 1000
    searchRadius = 5000
    
    # get the timestamp
    timeStamp = dt.datetime.strptime(str(combRad.iloc[0]['datetime']),'%Y-%m-%d %H:%M:%S')
    
    # generate the combined Total
    T = combineRadials(combRad, grid, searchRadius, bb[4], timeStamp)
 
 
    
    # add metadata related to bounding box
    #T = addBoundingBoxMetadata(T, lonMin, lonMax, latMin, latMax, gridResolution / 1000)
    
    # update is_combined attribute
    #T.is_combined = True
    
    print(T)
    T.plot()
    
    # add the Total object to the DataFrame
    combTot = pd.concat([combTot, pd.DataFrame([{'Total': T}])])
         
         
    # save file?
    
    return combTot

################################################################################


# testing
T = performRadialCombination(rDF, bb)

print(T)








#-------------------------------------------------------------------------------

# ok, I can now proceed to call functions





























    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

