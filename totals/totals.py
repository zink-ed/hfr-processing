import datetime as dt
from dateutil.relativedelta import relativedelta
import numpy as np
import xarray as xr
import netCDF4
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import re
import io
import os
from common import fileParser, addBoundingBoxMetadata
from collections import OrderedDict
from calc import evaluateGDOP, createLonLatGridFromBB
import fnmatch
import warnings
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from geopandas import GeoSeries
from statistics import mean


class Total(fileParser):
    
    """
    Totals Subclass

    this class should be used when loading CODAR (.tuv) total files
    """
    
    def __init__(self, fname='', replace_invalid=True, grid=gpd.GeoSeries(), empty_total=False):

        if not fname:
            empty_total = True
            replace_invalid = False
            
        super().__init__(fname)
        for key in self._tables.keys():
            table = self._tables[key]
            if 'LLUV' in table['TableType']:
                self.data = table['data']                
            elif 'src' in table['TableType']:
                self.diagnostics_source = table['data']
            elif 'CUR' in table['TableType']:
                self.cur_data = table['data']
                
        if 'SiteSource' in self.metadata.keys():
            if not self.is_wera:
                table_data = u''
                for ss in self.metadata['SiteSource']:
                    if '%%' in ss:
                        rep = {' comp': '_comp', ' Distance': '_Distance',' Ratio': '_Ratio',' (dB)': '_(dB)',' Width': '_Width', ' Resp': '_Resp', 'Value ': 'Value_','FOL ':'FOL_' }
                        rep = dict((re.escape(k), v) for k, v in rep.items())
                        pattern = re.compile('|'.join(rep.keys()))
                        ss_header = pattern.sub(lambda m: rep[re.escape(m.group(0))], ss).strip('%% SiteSource \n')
                    else:
                        ss = ss.replace('%SiteSource:', '').strip()
                        ss = ss.replace('Radial', '').strip()
                        table_data += '{}\n'.format(ss)
                # use pandas read_csv because it interprets the datatype for each column of the csv
                tdf = pd.read_csv(
                    io.StringIO(table_data),
                    sep=' ',
                    header=None,
                    names=ss_header.split(),
                    skipinitialspace=True
                )
                self.site_source = tdf
                
        # Evaluate GDOP for total files
        if hasattr(self, 'site_source'):
            # Get the site coordinates
            siteLon = self.site_source['Lon'].values.tolist()
            siteLat = self.site_source['Lat'].values.tolist()  
            # Create Geod object according to the Total CRS, if defined. Otherwise use WGS84 ellipsoid
            if self.metadata['GreatCircle']:
                g = Geod(ellps=self.metadata['GreatCircle'].split()[0].replace('"',''))                  
            else:
                g = Geod(ellps='WGS84')
                self.metadata['GreatCircle'] = '"WGS84"' + ' ' + str(g.a) + '  ' + str(1/g.f)
            self.data['GDOP'] = self.data.loc[:,['LOND','LATD']].apply(lambda x: evaluateGDOP(x,siteLon, siteLat, g),axis=1)                
        elif hasattr(self, 'data'):
            self.data['GDOP'] = np.nan
            
        # Evaluate the number of contributing radials (NRAD) for CODAR total files
        if hasattr(self, 'is_wera'):
            if not self.is_wera:
                if hasattr(self, 'data'):
                    self.data['NRAD'] = self.data.loc[:, self.data.columns.str.contains('S.*CN')].sum(axis=1)
            

        if replace_invalid:
            self.replace_invalid_values()
            
        if empty_total:
            self.empty_total()

        if not grid.empty:
            self.initialize_grid(grid)
            
    
    def empty_total(self):
        
        """
        create an empty Total object  by setting the geographical grid
        """

        self.file_path = ''
        self.file_name = ''
        self.full_file = ''
        self.metadata = ''
        self._iscorrupt = False
        self.time = []
        self.grid = GeoSeries()

        for key in self._tables.keys():
            table = self._tables[key]
            self._tables[key]['TableRows'] = '0'
            if 'LLUV' in table['TableType']:
                self.data.drop(self.data.index[:], inplace=True)
                self._tables[key]['data'] = self.data
            elif 'rads' in table['TableType']:
                self.diagnostics_radial.drop(self.diagnostics_radial.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_radial
            elif 'rcvr' in table['TableType']:
                self.diagnostics_hardware.drop(self.diagnostics_hardware.index[:], inplace=True)
                self._tables[key]['data'] = self.diagnostics_hardware
            elif 'RINF' in table['TableType']:
                self.range_information.drop(self.range_information.index[:], inplace=True)
                self._tables[key]['data'] = self.range_information
            elif 'CUR' in table['TableType']:
                self.cur_data.drop(self.cur_data.index[:], inplace=True)
                self._tables[key]['data'] = self.cur_data
                
        if not hasattr(self, 'data'):
            self.data = pd.DataFrame()
        
        if hasattr(self, 'site_source'):
            self.site_source.drop(self.site_source.index[:], inplace=True)
        else:
            self.site_source = pd.DataFrame()
            
    
    def initialize_grid(self, gridGS):
        
        """
        initialize the geogprahic grid for filling the LOND and LATD columns of the 
        Total object data DataFrame
        
        INPUT: gridGS = GeoPandas GeoSeries containing the longitude/latitude pairs of all
               the points in the grid
                
        OUTPUT: DataFrame with filled LOND and LATD columns
        """
        
        # initialize data DataFrame with column names
        self.data = pd.DataFrame(columns=['LOND', 'LATD', 'VELU', 'VELV', 'VELO', 'HEAD', 'UQAL', 'VQAL', 'CQAL', 'GDOP', 'NRAD'])
        
        # extract longitudes and latitude from grid GeoSeries and insert them into data DataFrame
        self.data['LOND'] = gridGS.x
        self.data['LATD'] = gridGS.y
        
        # add metadata about datum and CRS
        self.metadata = OrderedDict()
        self.metadata['GreatCircle'] = ''.join(gridGS.crs.ellipsoid.name.split()) + ' ' + str(gridGS.crs.ellipsoid.semi_major_metre) + '  ' + str(gridGS.crs.ellipsoid.inverse_flattening)
        
        
# plot - plot totals -----------------------------------------------------------
    def plot(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, shade=False, show=True, save=False, interpolated=False):
        
        """
        this function plots the current total velocity field (VELU and VELV components) on a 
        Cartesian grid
        
        The grid is defined either from the input values or from the Total object
        metadata. If no input is passed and no metadata related to the bounding box are present, the
        grid is defined from data content (LOND and LATD values).
        
        If 'shade' is False (default), a quiver plot with color and magnitude of the vectors proportional to
        current velocity is produced. If 'shade' is True, a quiver plot with uniform vetor lenghts is produced,
        superimposed to a pseudo-color map representing velocity magnitude.
        
        INPUT: lon_min = min lon value (decimal degrees)
               lon_max = max lon value (decimal degrees)
               lat_min = min lat value (decimal degrees)
               lat_max = max lat value (decimal degrees)
                 (if None it is taken from Total metadata)
               shade = boolean for enabling/disabling shade plot (default False)
               show = boolean for enabling/disabling plot visualization (default True)
            
        OUTPUT:

        """
        
        # initialize figure
        fig = plt.figure(figsize=(24, 16),tight_layout = {'pad': 0})
        #fig = plt.figure()
        
        # get the bounding box limits
        if not lon_min:
            lon_min = self.data.LOND.min() - 1
        if not lon_max:
            lon_max = self.data.LOND.max() + 1
        if not lat_min:
            lat_min = self.data.LATD.min() - 1 
        if not lat_max:
            lat_max = self.data.LATD.max() + 1   
                
        # evaluate lon and latof the center of the map
        lon_center = (lon_max + lon_min) / 2
        lat_center = (lat_max + lat_min) / 2
            
        # set the background map
        m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max, lon_0=lon_center, lat_0=lat_center, resolution = 'i', ellps='WGS84', projection='tmerc')        
        m.drawcoastlines()
        
        m.fillcontinents(color='#cc9955', lake_color='white')
        m.fillcontinents()
        m.drawparallels(np.arange(lat_min,lat_max))
        m.drawmeridians(np.arange(lon_min,lon_max))
        m.drawmapboundary(fill_color='white')
        
        m.bluemarble()
        
        # get station coordinates and codes
        siteLon = self.site_source['Lon'].to_numpy()
        siteLat = self.site_source['Lat'].to_numpy()
        siteCode = self.site_source['Name'].tolist()
        
        # compute the native map projection coordinates for the stations
        xS, yS = m(siteLon, siteLat)       
        
        print(siteCode)
        
        # plot radial stations
        m.plot(xS,yS,'rD')
        for label, xs, ys in zip(siteCode,xS,yS):
            plt.text(xs,ys,label,fontdict={'fontsize': 16, 'fontweight' : 'bold'})
        
        # plot velocity field
        if shade:
            '''
            self.to_xarray()
            
            # Create grid from longitude and latitude
            [longitudes, latitudes] = np.meshgrid(self.xdr['LONGITUDE'].data, self.xdr['LATITUDE'].data)
            
            # Compute the native map projection coordinates for the pseudo-color cells
            X, Y = m(longitudes, latitudes)

            # Create velocity variable in the shape of the grid
            V = abs(self.xdr['VELO'][0,0,:,:].to_numpy())
            # V = V[:-1,:-1]            
            
            # Make the pseudo-color plot
            warnings.simplefilter("ignore", category=UserWarning)
            c = m.pcolormesh(X, Y, V, shading='nearest', cmap=plt.cm.jet, vmin=0, vmax=1)
            '''
            from oceans.ocfis import uv2spdir, spdir2uv
            
            # Compute the native map projection coordinates for the vectors
            x, y = m(self.data.LOND, self.data.LATD)
            
            
            
            # Create the velocity component variables
            u = self.data.VELU / 100       # CODAR velocities are in cm/s
            v = self.data.VELV / 100       # CODAR velocities are in cm/s
            vel = abs(self.data.VELO) / 100
            
            angle, speed = uv2spdir(u.squeeze(), v.squeeze())
            u_norm, v_norm = spdir2uv(np.ones_like(speed), angle, deg=True)
            
            #color_clipped = np.clip(speed, 0, 1).squeeze()
            
            # Make the quiver plot
            m.quiver(x, y, u * 0.8 + u_norm * 0.2, v * 0.8 + v_norm * 0.2, vel, cmap=plt.cm.jet, width=0.001, headwidth=4, headlength=4, headaxislength=4)
            
            #warnings.simplefilter("default", category=UserWarning)
            
        else:
            # Compute the native map projection coordinates for the vectors
            x, y = m(self.data.LOND, self.data.LATD)
            
            # Create the velocity variables
            u = self.data.VELU / 100                # CODAR velocities are in cm/s
            v = self.data.VELV / 100                # CODAR velocities are in cm/s
            vel = abs(self.data.VELO) / 100         # CODAR velocities are in cm/s                
            
            # Make the quiver plot
            m.quiver(x, y, u, v, vel, cmap=plt.cm.jet, width=0.001, headwidth=4, headlength=4, headaxislength=4)
            
        # Add colorbar
        cbar = plt.colorbar()
        cbar.set_label('m/s',fontsize='x-large')
        
        # Add title
        plt.title(self.file_name + ' total velocity field', fontdict={'fontsize': 30, 'fontweight' : 'bold'})
                
        if show:
            plt.show()
        
        save_dir = '../total_images/'
        photo_name = self.file_name + '.png'
            
        if interpolated:
            print("INTERPOLATED")
            save_dir = save_dir + 'interpolated/'
            
        if save:
            print("SAVE")
            print(save_dir)
            print(photo_name)
            fig.savefig(save_dir + photo_name)
        
        return fig
        
        
    import numpy as np
    import datetime as dt
    import xarray as xr

    # def to_xarray_multidimensional
    def to_xarray(self, lon_min=None, lon_max=None, lat_min=None, lat_max=None, grid_res=None):
        
        """
        this function creates a dictionary of xarray DataArrays
          
        coordinate axes are set as (TIME, DEPTH, LATITUDE, LONGITUDE)
        
        INPUT: lon_min (decimal degrees)
               lat_max (decimal degrees)
               lat_min (decimal degrees)
               lat_max (decimal degrees)
               grid_res (meters)
        
            (otherwise info taken from metadata)
                
            
        OUTPUT: (generated dictionary is attached to Total object = xdr)

        """
        
        # initialize empty dictionary
        xdr = OrderedDict()
        
        '''    
        # longitude limits   
        if lon_min is None:
            lon_min = self.data.LOND.min()
        if lon_max is None:
            lon_max = self.data.LOND.max()
                
        # latitude limits   
        if lat_min is None:
            lat_min = self.data.LATD.min()
        if lat_max is None:
            lat_max = self.data.LATD.max()       
                
        # grid resolution
        if grid_res is None:                
            grid_res = float(3000)            
                
        # generate grid coordinates
        gridGS = createGrid(lon_min, lon_max, lat_min, lat_max, grid_res)
        '''   
                         
        # extract lons and lats from grid GeoSeries and insert them into numpy arrays
        lon_dim = np.unique(self.grid.x.to_numpy())
        lat_dim = np.unique(self.grid.y.to_numpy())
        
        # manage antimeridian crossing
        lon_dim = np.concatenate((lon_dim[lon_dim >= 0],lon_dim[lon_dim < 0]))         

        # create total grid from lons and lats
        [longitudes, latitudes] = np.meshgrid(lon_dim, lat_dim)

        # find grid indices from lon / lat grid (longitudes, latitudes)
        lat_map_idx = np.tile(np.nan, self.data['LATD'].shape)
        lon_map_idx = np.tile(np.nan, self.data['LOND'].shape)

        for i, line in enumerate(self.data['LATD']):
            lat_map_idx[i] = np.argmin(np.abs(lat_dim - self.data.LATD[i]))
            lon_map_idx[i] = np.argmin(np.abs(lon_dim - self.data.LOND[i]))
            
        # set X and Y coordinate mappings
        X_map_idx = lon_map_idx             # LONGITUDE is X axis
        Y_map_idx = lat_map_idx             # LATITUDE is Y axis
            
        # create dictionary containing variables from dataframe in the shape of total grid
        d = {key: np.tile(np.nan, longitudes.shape) for key in self.data.keys()}       
            
        # remap all variables
        for k, v in d.items():
            v[Y_map_idx.astype(int), X_map_idx.astype(int)] = self.data[k]
            d[k] = v

        # add extra dimensions for time (T) and depth (Z) - CF Standard: T, Z, Y, X -> T=axis0, Z=axis1
        d = {k: np.expand_dims(np.float32(v), axis=(0, 1)) for (k, v) in d.items()}            

        # drop LOND and LATD variables (they are set as coordinates of the DataSet)
        d.pop('LOND')
        d.pop('LATD')

        # scale velocities to be in m/s
        toMs = ['VELU', 'VELV', 'VELO', 'UQAL', 'VQAL','CQAL']
        for t in toMs:
            if t in d:
                d[t] = d[t] * 0.01

        # evaluate timestamp as number of days since 1950-01-01T00:00:00Z
        timeDelta = self.time - dt.datetime.strptime('1950-01-01T00:00:00Z','%Y-%m-%dT%H:%M:%SZ')
        ncTime = timeDelta.days + timeDelta.seconds / (60 * 60 * 24)
        
        # add all variables as xarray
        for k, v in d.items():
            xdr[k] = xr.DataArray(v, dims={'TIME': v.shape[0], 
                                           'DEPTH': v.shape[1], 
                                           'LATITUDE': v.shape[2], 
                                           'LONGITUDE': v.shape[3]}, 
                                     coords={'TIME': [ncTime], 
                                             'DEPTH': [0], 
                                             'LATITUDE': lat_dim, 
                                             'LONGITUDE': lon_dim})
        
        # add DataArray for coordinate variables
        xdr['TIME'] = xr.DataArray(ncTime,
                                   dims={'TIME': len(pd.date_range(self.time, periods=1))},
                                   coords={'TIME': [ncTime]})
        xdr['DEPTH'] = xr.DataArray(0,
                                    dims={'DEPTH': 1},
                                    coords={'DEPTH': [0]})
        xdr['LATITUDE'] = xr.DataArray(lat_dim,
                                       dims={'LATITUDE': lat_dim},
                                       coords={'LATITUDE': lat_dim})
        xdr['LONGITUDE'] = xr.DataArray(lon_dim,
                                        dims={'LONGITUDE': lon_dim},
                                        coords={'LONGITUDE': lon_dim})  
        
        # attach the dictionary to the Total object
        self.xdr = xdr
        
        return
    
    
    def check_ehn_mandatory_variables(self):
        """
        This function checks if the Total object contains all the mandatory data variables
        (i.e. not coordinate variables) required by the European standard data model developed in the framework of the 
        EuroGOOS HFR Task Team.
        Missing variables are appended to the DataFrame containing data, filled with NaNs.
        
        INPUT:            
            
        OUTPUT:
        """
        # Set mandatory variables based on the HFR manufacturer
        chkVars = ['VELU', 'VELV', 'UQAL', 'VQAL', 'CQAL', 'GDOP']
            
        # Check variables and add missing ones
        for vv in chkVars:
            if vv not in self.data.columns:
                self.data[vv] = np.nan
                
        return
    
    
    def apply_ehn_datamodel(self, network_data, station_data, version):
        """
        This function applies the European standard data model developed in the
        framework of the EuroGOOS HFR Task Team to the Total object.
        The Total object content is stored into an xarray Dataset built from the
        xarray DataArrays created by the Total method to_xarray_multidimensional.
        Variable data types and data packing information are collected from
        "Data_Models/EHN/Totals/Total_Data_Packing.json" file.
        Variable attribute schema is collected from 
        "Data_Models/EHN/Totals/Total_Variables.json" file.
        Global attribute schema is collected from 
        "Data_Models/EHN/Global_Attributes.json" file.
        Global attributes are created starting from Total object metadata and from 
        DataFrames containing the information about HFR network and radial stations
        read from the EU HFR NODE database.
        The generated xarray Dataset is attached to the Total object, named as xds.
        
        INPUT:
            
            
        OUTPUT:
        """
        # Set the netCDF format
        ncFormat = 'NETCDF4_CLASSIC'
        
        # Expand Total object variables along the coordinate axes
        if self.is_combined:
            self.to_xarray_multidimensional()
        else:
            # Get bounding box limits and grid resolution from database
            lonMin = network_data.iloc[0]['geospatial_lon_min']
            lonMax = network_data.iloc[0]['geospatial_lon_max']
            latMin = network_data.iloc[0]['geospatial_lat_min']
            latMax = network_data.iloc[0]['geospatial_lat_max']
            gridRes = network_data.iloc[0]['grid_resolution']*1000
            self.to_xarray_multidimensional(lonMin,lonMax,latMin,latMax,gridRes)
        
        # Set auxiliary coordinate sizes
        maxsiteSize = 150
        refmaxSize = 50
        maxinstSize = 50
        
        # Get data packing information per variable
        f = open('Data_Models/EHN/Totals/Total_Data_Packing.json')
        dataPacking = json.loads(f.read())
        f.close()
        
        # Get variable attributes
        f = open('Data_Models/EHN/Totals/Total_Variables.json')
        totVariables = json.loads(f.read())
        f.close()
        
        # Get global attributes
        f = open('Data_Models/EHN/Global_Attributes.json')
        globalAttributes = json.loads(f.read())
        f.close()
        
        # Rename velocity related and quality related variables
        self.xdr['EWCT'] = self.xdr.pop('VELU')
        self.xdr['NSCT'] = self.xdr.pop('VELV')
        if 'UQAL' in self.xdr:
            self.xdr['EWCS'] = self.xdr.pop('UQAL')
        if 'VQAL' in self.xdr:
            self.xdr['NSCS'] = self.xdr.pop('VQAL')
        if 'CQAL' in self.xdr:
            self.xdr['CCOV'] = self.xdr.pop('CQAL')
        
        # Drop unnecessary DataArrays from the DataSet
        toDrop = ['VFLG', 'XDST', 'YDST', 'RNGE', 'BEAR','NRAD', 'VELO', 'HEAD','index']
        for t in toDrop:
            if t in self.xdr:
                self.xdr.pop(t)
        toDrop = list(self.xdr.keys())
        for t in toDrop:
            if fnmatch.fnmatch(t,'S*CN'):
                self.xdr.pop(t)
        toDrop = []
        for vv in self.xdr:
            if vv not in totVariables.keys():
                toDrop.append(vv)
        for rv in toDrop:
            self.xdr.pop(rv)            
            
        # Add coordinate reference system to the dictionary
        self.xdr['crs'] = xr.DataArray(int(0), )       
        
        # Add antenna related variables to the dictionary
        # Number of antennas        
        contributingSiteNrx = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())]['number_of_receive_antennas'].to_numpy()
        nRX = np.asfarray(contributingSiteNrx)
        nRX = np.pad(nRX, (0, maxsiteSize - len(nRX)), 'constant',constant_values=(np.nan,np.nan))
        contributingSiteNtx = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())]['number_of_transmit_antennas'].to_numpy()
        nTX = np.asfarray(contributingSiteNtx)
        nTX = np.pad(nTX, (0, maxsiteSize - len(nTX)), 'constant',constant_values=(np.nan,np.nan))
        self.xdr['NARX'] = xr.DataArray([nRX], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['NATX'] = xr.DataArray([nTX], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        
        # Longitude and latitude of antennas
        contributingSiteLat = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())]['site_lat'].to_numpy()
        siteLat = np.pad(contributingSiteLat, (0, maxsiteSize - len(contributingSiteLat)), 'constant',constant_values=(np.nan,np.nan))
        contributingSiteLon = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())]['site_lon'].to_numpy()
        siteLon = np.pad(contributingSiteLon, (0, maxsiteSize - len(contributingSiteLon)), 'constant',constant_values=(np.nan,np.nan))
        self.xdr['SLTR'] = xr.DataArray([siteLat], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNR'] = xr.DataArray([siteLon], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLTT'] = xr.DataArray([siteLat], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SLNT'] = xr.DataArray([siteLon], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        
        # Codes of antennas
        contributingSiteCodeList = station_data.loc[station_data['station_id'].isin(self.site_source.Name.tolist())]['station_id'].tolist()
        antCode = np.array([site.encode() for site in contributingSiteCodeList])
        antCode = np.pad(antCode, (0, maxsiteSize - len(contributingSiteCodeList)), 'constant',constant_values=('',''))
        self.xdr['SCDR'] = xr.DataArray(np.array([antCode]), dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDR'].encoding['char_dim_name'] = 'STRING' + str(len(station_data['station_id'].to_numpy()[0]))
        self.xdr['SCDT'] = xr.DataArray(np.array([antCode]), dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXSITE': maxsiteSize})
        self.xdr['SCDT'].encoding['char_dim_name'] = 'STRING' + str(len(station_data['station_id'].to_numpy()[0]))
                
        # Add SDN namespace variables to the dictionary
        siteCode = ('%s' % network_data.iloc[0]['network_id']).encode()
        self.xdr['SDN_CRUISE'] = xr.DataArray([siteCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_CRUISE'].encoding['char_dim_name'] = 'STRING' + str(len(siteCode))
        platformCode = ('%s' % network_data.iloc[0]['network_id'] + '-Total').encode()
        self.xdr['SDN_STATION'] = xr.DataArray([platformCode], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_STATION'].encoding['char_dim_name'] = 'STRING' + str(len(platformCode))
        ID = ('%s' % platformCode.decode() + '_' + self.time.strftime('%Y-%m-%dT%H:%M:%SZ')).encode()
        self.xdr['SDN_LOCAL_CDI_ID'] = xr.DataArray([ID], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_LOCAL_CDI_ID'].encoding['char_dim_name'] = 'STRING' + str(len(ID))
        sdnEDMO = np.asfarray(pd.concat([network_data['EDMO_code'],station_data['EDMO_code']]).unique())
        sdnEDMO = np.pad(sdnEDMO, (0, maxinstSize - len(sdnEDMO)), 'constant',constant_values=(np.nan,np.nan))
        self.xdr['SDN_EDMO_CODE'] = xr.DataArray([sdnEDMO], dims={'TIME': len(pd.date_range(self.time, periods=1)), 'MAXINST': maxinstSize})
        sdnRef = ('%s' % network_data.iloc[0]['metadata_page']).encode()
        self.xdr['SDN_REFERENCES'] = xr.DataArray([sdnRef], dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['SDN_REFERENCES'].encoding['char_dim_name'] = 'STRING' + str(len(sdnRef))
        sdnXlink = ('%s' % '<sdn_reference xlink:href=\"' + sdnRef.decode() + '\" xlink:role=\"\" xlink:type=\"URL\"/>').encode()
        self.xdr['SDN_XLINK'] = xr.DataArray(np.array([[sdnXlink]]), dims={'TIME': len(pd.date_range(self.time, periods=1)), 'REFMAX': refmaxSize})
        self.xdr['SDN_XLINK'].encoding['char_dim_name'] = 'STRING' + str(len(sdnXlink))
        
        # Add spatial and temporal coordinate QC variables (set to good data due to the nature of HFR system)
        self.xdr['TIME_QC'] = xr.DataArray([1],dims={'TIME': len(pd.date_range(self.time, periods=1))})
        self.xdr['POSITION_QC'] = self.xdr['QCflag'] * 0 + 1
        self.xdr['DEPTH_QC'] = xr.DataArray([1],dims={'TIME': len(pd.date_range(self.time, periods=1))})
            
        # Create DataSet from DataArrays
        self.xds = xr.Dataset(self.xdr)
        
        # Add data variable attributes to the DataSet
        for vv in self.xds:
            self.xds[vv].attrs = totVariables[vv]
            
        # Update QC variable attribute "comment" for inserting test thresholds and attribute "flag_values" for assigning the right data type
        for qcv in self.metadata['QCTest']:
            if qcv in self.xds:
                self.xds[qcv].attrs['comment'] = self.xds[qcv].attrs['comment'] + ' ' + self.metadata['QCTest'][qcv]
                self.xds[qcv].attrs['flag_values'] = list(np.int_(self.xds[qcv].attrs['flag_values']).astype(dataPacking[qcv]['dtype']))
        for qcv in ['TIME_QC', 'POSITION_QC', 'DEPTH_QC']:
            if qcv in self.xds:
                self.xds[qcv].attrs['flag_values'] = list(np.int_(self.xds[qcv].attrs['flag_values']).astype(dataPacking[qcv]['dtype']))
                
        # Add coordinate variable attributes to the DataSet
        for cc in self.xds.coords:
            self.xds[cc].attrs = totVariables[cc]
            
        # Evaluate measurement maximum depth
        vertMax = 3e8 / (8*np.pi * station_data['transmit_central_frequency'].to_numpy().min()*1e6)
        
        # Evaluate time coverage start, end, resolution and duration
        timeCoverageStart = self.time - relativedelta(minutes=network_data.iloc[0]['temporal_resolution']/2)
        timeCoverageEnd = self.time + relativedelta(minutes=network_data.iloc[0]['temporal_resolution']/2)
        timeResRD = relativedelta(minutes=network_data.iloc[0]['temporal_resolution'])
        timeCoverageResolution = 'PT'
        if timeResRD.hours !=0:
            timeCoverageResolution += str(int(timeResRD.hours)) + 'H'
        if timeResRD.minutes !=0:
            timeCoverageResolution += str(int(timeResRD.minutes)) + 'M'
        if timeResRD.seconds !=0:
            timeCoverageResolution += str(int(timeResRD.seconds)) + 'S'   
            
        # Fill global attributes
        globalAttributes['site_code'] = siteCode.decode()
        globalAttributes['platform_code'] = platformCode.decode()
        globalAttributes.pop('oceanops_ref')
        globalAttributes.pop('wmo_platform_code')
        globalAttributes.pop('wigos_id')
        globalAttributes['doa_estimation_method'] = ', '.join(station_data[["station_id", "DoA_estimation_method"]].apply(": ".join, axis=1))
        globalAttributes['calibration_type'] = ', '.join(station_data[["station_id", "calibration_type"]].apply(": ".join, axis=1))
        if 'HFR-US' in network_data.iloc[0]['network_id']:
            station_data['last_calibration_date'] = 'N/A'
            globalAttributes['last_calibration_date'] = ', '.join(pd.concat([station_data['station_id'],station_data['last_calibration_date']],axis=1)[["station_id", "last_calibration_date"]].apply(": ".join, axis=1))
        else:
            globalAttributes['last_calibration_date'] = ', '.join(pd.concat([station_data['station_id'],station_data['last_calibration_date'].apply(lambda x: x.strftime('%Y-%m-%dT%H:%M:%SZ'))],axis=1)[["station_id", "last_calibration_date"]].apply(": ".join, axis=1))
            globalAttributes['last_calibration_date'] = globalAttributes['last_calibration_date'].replace('1-01-01T00:00:00Z', 'N/A')
        globalAttributes['calibration_link'] = ', '.join(station_data[["station_id", "calibration_link"]].apply(": ".join, axis=1))
        # globalAttributes['title'] = network_data.iloc[0]['title']
        globalAttributes['title'] = 'Near Real Time Surface Ocean Total Velocity by ' + globalAttributes['platform_code']
        globalAttributes['summary'] = network_data.iloc[0]['summary']
        globalAttributes['institution'] = ', '.join(pd.concat([network_data['institution_name'],station_data['institution_name']]).unique().tolist())
        globalAttributes['institution_edmo_code'] = ', '.join([str(x) for x in pd.concat([network_data['EDMO_code'],station_data['EDMO_code']]).unique().tolist()])
        globalAttributes['institution_references'] = ', '.join(pd.concat([network_data['institution_website'],station_data['institution_website']]).unique().tolist())
        globalAttributes['id'] = ID.decode()
        globalAttributes['project'] = network_data.iloc[0]['project']
        globalAttributes['comment'] = network_data.iloc[0]['comment']
        globalAttributes['network'] = network_data.iloc[0]['network_name']
        globalAttributes['data_type'] = globalAttributes['data_type'].replace('current data', 'total current data')
        globalAttributes['geospatial_lat_min'] = str(network_data.iloc[0]['geospatial_lat_min'])
        globalAttributes['geospatial_lat_max'] = str(network_data.iloc[0]['geospatial_lat_max'])
        globalAttributes['geospatial_lat_resolution'] = str(network_data.iloc[0]['grid_resolution'])
        globalAttributes['geospatial_lon_min'] = str(network_data.iloc[0]['geospatial_lon_min'])
        globalAttributes['geospatial_lon_max'] = str(network_data.iloc[0]['geospatial_lon_max'])
        globalAttributes['geospatial_lon_resolution'] = str(network_data.iloc[0]['grid_resolution'])        
        globalAttributes['geospatial_vertical_max'] = str(vertMax)
        globalAttributes['geospatial_vertical_resolution'] = str(vertMax)        
        globalAttributes['time_coverage_start'] = timeCoverageStart.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['time_coverage_end'] = timeCoverageEnd.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['time_coverage_resolution'] = timeCoverageResolution
        globalAttributes['time_coverage_duration'] = timeCoverageResolution
        globalAttributes['area'] = network_data.iloc[0]['area']
        globalAttributes['format_version'] = version
        globalAttributes['netcdf_format'] = ncFormat
        globalAttributes['citation'] += network_data.iloc[0]['citation_statement']
        globalAttributes['license'] = network_data.iloc[0]['license']
        globalAttributes['acknowledgment'] = network_data.iloc[0]['acknowledgment']
        globalAttributes['processing_level'] = '3B'
        globalAttributes['contributor_name'] = network_data.iloc[0]['contributor_name']
        globalAttributes['contributor_role'] = network_data.iloc[0]['contributor_role']
        globalAttributes['contributor_email'] = network_data.iloc[0]['contributor_email']
        globalAttributes['manufacturer'] = ', '.join(station_data[["station_id", "manufacturer"]].apply(": ".join, axis=1))
        globalAttributes['sensor_model'] = ', '.join(station_data[["station_id", "manufacturer"]].apply(": ".join, axis=1))
        globalAttributes['software_version'] = version
        
        creationDate = dt.datetime.utcnow()
        globalAttributes['metadata_date_stamp'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_created'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['date_modified'] = creationDate.strftime('%Y-%m-%dT%H:%M:%SZ')
        globalAttributes['history'] = 'Data collected at ' + self.time.strftime('%Y-%m-%dT%H:%M:%SZ') + '. netCDF file created at ' \
                                    + creationDate.strftime('%Y-%m-%dT%H:%M:%SZ') + ' by the European HFR Node.'        
        
        # Add global attributes to the DataSet
        self.xds.attrs = globalAttributes
            
        # Encode data types, data packing and _FillValue for the data variables of the DataSet
        for vv in self.xds:
            if vv in dataPacking:
                if 'dtype' in dataPacking[vv]:
                    self.xds[vv].encoding['dtype'] = dataPacking[vv]['dtype']
                if 'scale_factor' in dataPacking[vv]:
                    self.xds[vv].encoding['scale_factor'] = dataPacking[vv]['scale_factor']                
                if 'add_offset' in dataPacking[vv]:
                    self.xds[vv].encoding['add_offset'] = dataPacking[vv]['add_offset']
                if 'fill_value' in dataPacking[vv]:
                    self.xds[vv].encoding['_FillValue'] = netCDF4.default_fillvals[np.dtype(dataPacking[vv]['dtype']).kind + str(np.dtype(dataPacking[vv]['dtype']).itemsize)]
                else:
                    self.xds[vv].encoding['_FillValue'] = None
                    
        # Update valid_min and valid_max variable attributes according to data packing
        for vv in self.xds:
            if 'valid_min' in totVariables[vv]:
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_min'] = np.float_(((totVariables[vv]['valid_min'] - dataPacking[vv]['add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_min'] = np.float_(totVariables[vv]['valid_min']).astype(dataPacking[vv]['dtype'])
            if 'valid_max' in totVariables[vv]:             
                if ('scale_factor' in dataPacking[vv]) and ('add_offset' in dataPacking[vv]):
                    self.xds[vv].attrs['valid_max'] = np.float_(((totVariables[vv]['valid_max'] - dataPacking[vv]['add_offset']) / dataPacking[vv]['scale_factor'])).astype(dataPacking[vv]['dtype'])
                else:
                    self.xds[vv].attrs['valid_max'] = np.float_(totVariables[vv]['valid_max']).astype(dataPacking[vv]['dtype'])
            
        # Encode data types and avoid data packing, valid_min, valid_max and _FillValue for the coordinate variables of the DataSet
        for cc in self.xds.coords:
            if cc in dataPacking:
                if 'dtype' in dataPacking[cc]:
                    self.xds[cc].encoding['dtype'] = dataPacking[cc]['dtype']
                if 'valid_min' in totVariables[cc]:
                    del self.xds[cc].attrs['valid_min']
                if 'valid_max' in totVariables[cc]:
                    del self.xds[cc].attrs['valid_max']
                self.xds[cc].encoding['_FillValue'] = None
               
        return

    
    def initialize_qc(self):
        
        """
        initialize dictionary entry for QC metadata
        """
        
        # initialize dictionary entry for QC metadta
        self.metadata['QCTest'] = {}
        
        
        
        
    # qc_ehn_data_density_threshold - QC301    
    def qc_qartod_data_density_threshold(self, minContrRad=2):
        
        """
        this test labels total velocity vectors with a number of contributing radial velocities bigger 
        than the minimum number defined for normal operations with a “good data” flag
        
        INPUTS: minContrRad = min number of contributing radial velocities for normal operations                     
        """
        
        # set the test name
        testName = 'QC301'
        
        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:,testName] = 1
    
        # set bad flag for velocities not passing the test
        if 'NRAD' in self.data.columns:
            self.data.loc[(self.data['NRAD'] < minContrRad), testName] = 4
            print("DATA DENSITY THRESHOLD")
    
        self.metadata['QCTest'][testName] = 'Data Density Threshold QC Test - Test applies to each vector. ' \
            + 'Threshold=[' + f'minimum number of contributing radial velocities={minContrRad}]'
        
    
    
        # qc_ehn_threshold - QC302    
    def qc_qartod_gdop_threshold(self, maxGDOP=2):
        
        """
        this test labels total velocity vectors whose GDOP is smaller than a maximum GDOP threshold 
        with a “good data” flag
        
        INPUTS: maxGDOP = maximum allowed GDOP for normal operations                     
        """
        
        # set the test name
        testName = 'QC302'
        
        # add new column to the DataFrame, set every row as passing the test (flag = 1)
        self.data.loc[:,testName] = 1
    
        # set bad flag for velocities not passing the test
        self.data.loc[(self.data['GDOP'] > maxGDOP), testName] = 4
    
        self.metadata['QCTest'][testName] = 'GDOP Threshold QC Test - Test applies to each vector. ' \
            + 'Threshold=[' + f'GDOP threshold={maxGDOP}]'
        
     
     
        
    # qc_ehn_maximum_velocity - QC303    
    def qc_qartod_maximum_velocity(self, maxSpeed=250):
        
        """
        this test labels total velocity vectors whose module is smaller than a maximum velocity threshold 
        with a “good data” flag
        
        INPUTS: maxSpeed = maximum velocity in m/s for normal operations                     
        """
        
        # set the test name
        testName = "QC303"
        
        # make sure VELO is float
        self.data["VELO"] = self.data["VELO"].astype(float)
        
        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        self.data.loc[:,testName] = 1
    
        # set bad flag for velocities not passing the test
        self.data.loc[(self.data['VELO'].abs() > maxSpeed), testName] = 4
    
        self.metadata['QCTest'][testName] = 'Velocity Threshold QC Test - Test applies to each vector. ' \
            + 'Threshold=[' + f'maximum velocity={maxSpeed} (cm/s)]'

       
    # QC304
    def qc_qartod_spatial_median(self, dx, dy, limit = 2, diff=25):
        
        testName = 'QC304'
        
        self.data.loc[:,testName] = 1
        
        #print(dx)
        #print(dy)
        
        velo = self.data['VELO'].to_numpy()
        #print(velo.shape)
        diffcol = self.data['VELO'].to_numpy() + np.nan
        
        velo = velo.reshape((dx, dy))
        #print(velo.shape)
        diffcol = diffcol.reshape((dx, dy))
        
        
        padded = np.pad(velo, pad_width=((limit, limit), (limit, limit)), mode='constant', constant_values=np.nan)
        diffcol = np.pad(velo, pad_width=((limit, limit), (limit, limit)), mode='constant', constant_values=np.nan)
        #print(padded.shape)

        for row in range(limit, dx + limit):
            for col in range(limit, dy + limit):
                temp = padded[row - limit: row + limit + 1, col - limit: col + limit + 1]
                median = np.nanmedian(temp)
                diffcol[row][col] = padded[row][col] - median
        
        diffcol = diffcol[limit:-limit, limit:-limit]
        diffcol = diffcol.flatten()
        
        boolean = np.abs(diffcol) > diff
        
        
        self.data[testName] = self.data[testName].where(~boolean, other=4)
        
        print(self.data[self.data[testName] == 4])
        self.metadata['QCTest'][testName] = 'Spatial Median QC Test - Test applies to each vector. ' \
            + 'Thresholds=[' + f'current_difference={str(diff)} (cm/s)]'
        
    
    # qc_ehn_temporal_derivative - QC206    
    def qc_qartod_temporal_derivative(self, t0, tempDerThr=100):
        
        """
        this test compares the velocity of each total vector with the velocity of the total vector 
        measured in the previous timestamp at the same location.
        Each vector for which the velocity difference is smaller than the specified threshold for normal 
        operations (tempDerThr), is labeled with a "good data" flag.
        
        INPUTS: t0 = Total object of the previous timestamp
                tempDerThr = velocity difference threshold in m/s for normal operations
        """
        
        # Set the test name
        testName = 'QC206'
        
        # Check if the previous timestamp total file exists
        if not t0 is None:
            # Merge the data DataFrame of the two Totals and evaluate velocity differences at each location
            mergedDF = self.data.merge(t0.data, on=['LOND', 'LATD'], how='left', suffixes=(None, '_x'), indicator='Exist')
            velDiff = (mergedDF['VELO'] - mergedDF['VELO_x']).abs()

            # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
            self.data.loc[:,testName] = 1

            # Set rows of the DataFrame for QC data as not evaluated (flag = 0) for locations existing in the current total but not in the previous one
            self.data.loc[mergedDF['Exist'] == 'left_only', testName] = 0

            # Set bad flag for vectors not passing the test
            self.data.loc[(velDiff > tempDerThr), testName] = 4         # velocity in cm/s (LLUV)

        else:
            # Add new column to the DataFrame for QC data by setting every row as not evaluated (flag = 0)
            self.data.loc[:,testName] = 0
        
        self.metadata['QCTest'][testName] = 'Temporal Derivative QC Test - Test applies to each vector. ' \
            + 'Threshold=[' + f'velocity difference threshold={str(tempDerThr)} (cm/s)]'
 
    
    # qc_ehn_overall_qc_flag    
    def qc_qartod_primary_flag(self, erase=False):
        
        """
        this QC test labels total velocity vectors with a ‘good_data” flag if all QC tests are passed
        
        INPUTS: erase = if True, then it will remove data with a bad flag
            
        """
        
        # Set the test name
        testName = 'PRIM'
        
        # Add new column to the DataFrame for QC data by setting every row as not passing the test (flag = 4)
        self.data[testName] = 4
        
        # Set good flags for vectors passing all QC tests
        self.data.loc[self.data.loc[:, self.data.columns.str.contains('QC')].eq(1).all(axis=1), testName] = 1

        self.metadata['QCTest'][testName] = 'Overall QC Flag - Test applies to each vector. Test checks if all QC tests are passed.'
        
        #print(self.data)
        
        if erase:
            indexes = self.data[self.data['PRIM'] == 4].index
            self.data.drop(indexes, inplace=True)
            self.data.reset_index(level=None, drop=False, inplace=True)
            print("ERASED")
             

    def file_type(self):
        """
        Return a string representing the type of file this is.
        """
        return 'totals'
        
        
        
