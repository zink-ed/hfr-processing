
# Created Sep 2024 by Cathleen Qiao

# This file is for interpolating total vectors.

# Look at these repositories for more information:
# https://github.com/LorenzoCorgnati/EU_HFR_NODE_pyHFR/tree/main
# https://github.com/rucool/hfradarpy/tree/master

# This file is unfinished. 


import numpy as np
from statistics import mean

# check for row or column interpolation (have neighbors)   
def check_next(m, r, c, rl, cl):
    if (not m[r][c-1] and not m[r][c+1]):
        rl.append(r)
        rl.append(r)
        cl.append(c-1)
        cl.append(c+1)
    if (not m[r-1][c] and not m[r+1][c]):
        rl.append(r-1)
        rl.append(r+1)
        cl.append(c)
        cl.append(c)
    if (len(rl) > 0):
        return True
    return False
    
    
def calculate_mean(rl, cl, m):
    t = []
    for r, c in zip(rl, cl):
        t.append(m[r][c])
    result = mean(t)
    
    return result


# checking if grid is there
def check_grid(m, r, c, a):
    if (m[r + a[0]][c - a[2]] or 
        m[r - a[1]][c - a[2]] or
        m[r + a[0]][c + a[3]] or
        m[r - a[1]][c + a[3]]):
            return False
    return True


# possible grids
def get_grid(m, r, c, a):
    
    # check [1, 1, 1, 1]
    if check_grid(m, r, c, a):
        return True

    # check when one entry is 2
    for e in a:
        e = 2
        if check_grid(m, r, c, a):
            return True
        e = 1
    
    # when two entries are 2
    for i in range(4):
        a[i] = 2
        for j in range(i + 1, 4):
            a[j] = 2
            if check_grid(m, r, c, a):
                return True
            a[j] = 1
        a[i] = 1     
    
    return False
        
        
# bilinear interpolation formula
def bilinear(m, r, c, a):

    b00 = m[r + a[0]][c - a[2]]
    b01 = m[r - a[1]][c - a[2]]
    b10 = m[r + a[0]][c + a[3]]
    b11 = m[r - a[1]][c + a[3]]
    
    b1 = (b00 * a[0] + b01 * a[1]) / (a[0] + a[1])
    b2 = (b10 * a[0] + b11 * a[1]) / (a[0] + a[1])
    
    m[r][c] = (b1 * a[2] + b2 * a[3]) / (a[2] + a[3])
    


def applyFlag(T):
   
        # set the test name
        testName = 'INTP'
        
        # Add new column to the DataFrame for QC data by setting every row as passing the test (flag = 1)
        T.data.loc[:,testName] = 0
    
        # set flag for points with data
        T.data.loc[(T.data['VELU'].isnull() == False), testName] = 2
        T.data.loc[(T.data['HEAD'].isnull() == False), testName] = 1


 
def interpolation(T, l, w, d):
    
    #print(T.data)
    #print(T.data.dropna(subset=['VELU']))
    
    #lon_orig = T.data.LOND.to_numpy()
    #lat_orig = T.data.LATD.to_numpy()
    u_orig = T.data.VELU.to_numpy()
    v_orig = T.data.VELV.to_numpy()
    velo_orig = T.data.VELO.to_numpy()
    #print(u_orig.shape)
    
    u_orig = u_orig.reshape(l, w)
    v_orig = v_orig.reshape(l, w)
    velo_orig = velo_orig.reshape(l, w)
    original = np.isnan(u_orig)
    #print(u_orig.shape)
    #print(u_orig[~np.isnan(u_orig)].shape)
    
    for r in range(d, l - d):
        
        for c in range(d, w - d):
            
            a = [1, 1, 1, 1]
            
            # if lon + lat has been interpolated
            if original[r][c] == True:
                
                rl = []
                cl = []
            
                # for points that are close enough
                if (check_next(original, r, c, rl, cl)):
                    u = calculate_mean(rl, cl, u_orig)
                    #print(u_orig[r][c])
                    u_orig[r][c] = u
                    #print(u_orig[r][c])
                    v = calculate_mean(rl, cl, v_orig)
                    v_orig[r][c] = v
                    velo = calculate_mean(rl, cl, velo_orig)
                    velo_orig[r][c] = velo
                    #print("FIRST")
            
                # add points if grid has been found
                elif get_grid(original, r, c, a):
                    bilinear(u_orig, r, c, a)
                    bilinear(v_orig, r, c, a)
                    bilinear(velo_orig, r, c, a)
                    #print("SECOND")
    
    #print(u_orig[~np.isnan(u_orig)].shape)
    
    #print(u_orig.flatten().tolist())
    #T.data['VELU'] = u_orig.flatten().tolist()
    #T.data.VELV = v_orig.flatten().tolist()
    #original = original.flatten().tolist()
    
    #print(T.data.dropna(subset=['VELU']))
    #T.data.VELO = velo_orig.flatten()
    
    applyFlag(T)
    
    return T
    
    

