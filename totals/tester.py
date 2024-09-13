
# Created Sep 2024 by Cathleen Qiao

# This file is for testing functions by reading a pickle file. This is useful
# for testing as the total does not need to be recreated (the pickle file
# acts as a breakpoint). I use this file to test my interpolation functions
# and plotting functions.


from totals import Total
import total_interpolation as ti
import pickle


def applyQC(T, dx, dy):

    # check if Total object contains data
    if T.data.empty:
        print('total file is empty: no QC test applied')
        return T
        
    # initialize QC metadata
    T.initialize_qc()
    
    # DDNS
    T.qc_qartod_data_density_threshold()
    
    # CSPD
    T.qc_qartod_maximum_velocity()

    # Temporal Gradient
    '''
    prevHourTime = T.time-dt.timedelta(minutes=networkData.iloc[0]['temporal_resolution'])
    prevHourFileName = buildEHNtotalFilename(networkData.iloc[0]['network_id'],prevHourTime,'.ttl')
    prevHourTotFile = prevHourFolderPath + prevHourFileName     # previous hour total file
    if os.path.exists(prevHourTotFile):
        with open(prevHourTotFile, 'rb') as ttlFile:
            t0 = pickle.load(ttlFile)
    '''
    
    # GDOP
    T.qc_qartod_gdop_threshold()
    
    # Spatial Median
    T.qc_qartod_spatial_median(dx, dy)

    # Overall QC
    T.qc_qartod_primary_flag(True)
    
    return T
    
time = '2024_04_20_0000'
pkl = '../total-data/' + 'total_' + time + '.pkl'

with open(pkl, 'rb') as file:

    T = pickle.load(file)
    print(T.data)

#print(T.data.dropna(subset=['VELU']))
#print(T.data.dropna(subset=['HEAD']))
  
#T = ti.interpolation(T, 153, 76, 2)
#T.data = T.data.dropna(subset=['VELU']) 
#print(T.data)
    
#T = applyQC(T, 153, 76)
#print(T.data)

#T.data = T.data.dropna(subset=['VELU'])
#print(T.data)


    
T.plot_cartopy(show=False, save=True, interpolated=False)







