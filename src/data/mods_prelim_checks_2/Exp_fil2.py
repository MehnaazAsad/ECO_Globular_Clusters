# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 17:36:35 2017

@author: asadm2
"""
### DESCRIPTION
#This script carries out an exposure time check and an "at least 2 filters" 
#check 

import pandas as pd
import numpy as np

def expfil2(objname,revtxt):

    """
    This function carries out an exposure time check and also checks if 
    the images that passed the first check are taken using at least 2 filters
    
    Args: 
        objname: ECOID of the galaxy
        revtxt: objname_rev.txt file that Obj_in_Img.py returns
    
    Returns:
        goodObj.txt and badObj.txt files depending on which ECOIDs passed
        both checks and which ones didn't
    """
    path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
    'data/raw/'
    path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
    'data/interim/'
    
    ECO = path_to_raw + 'Available_HST_Data_ECO.txt'
    ECO = pd.read_csv(ECO, delimiter='\s+', header=None, \
                      names=['ECOID', 'HSTOBJ', 'RA', 'DEC', 'exptime', \
                             'camera', 'filename'])
    
    ECO['exptime'] = pd.to_numeric(ECO['exptime'],errors='coerce')
    ECO['filename'] = ECO['filename'].astype('str')
    
    files_arr = ECO['filename'].values
    n_files = len(files_arr)
    
    wfc3_ir = ['f110w','f125w','f160w']   
    wfc3_uvis = ['f606w','f600lp']
    filters = [[] for x in range(n_files)]
    for i in range(len(ECO['filename'])):
        str_split = ECO['filename'][i].split(';')[1].split('_')
        filter_i = str_split[3]+'_'+str_split[4]+'_'+str_split[5]
    
        if 'ACS' in filter_i: #acs_wfc
            filter_i = filter_i.lower()
            
        elif 'd634' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
            
        elif 'm51' in filter_i: #acs-wfc fixed to acs_wfc       
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'tile' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'c_v' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'ngc' in filter_i: #acs fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_wfc_'+str_split[2]
            
        elif '131009' in filter_i: #wfc3
            str_split = filter_i.split('_')
            if str_split[2] == 'f438w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            elif str_split[2] == 'f775w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
                                                 
        elif 'par' in filter_i:#and any(str in filter_i for str in wfc3_ir):
            str_split = filter_i.split('_')
            if str_split[2] == 'f606w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            elif str_split[2] == 'f125w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f160w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f110w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f600lp':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            
        elif 'w_wf' in filter_i: #wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[1]
        
        elif 'lp_wf' in filter_i: #wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[1]
            
        elif 'n4496' in filter_i: #all wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]            
            
        elif 'n5194' in filter_i: #all wfpc2 
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
           
        elif 'u6614' in filter_i: #all wfpc2   
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            
        filters[i] = filter_i
        
    filters = np.asarray(filters)
    
    filters_unique = np.unique(filters)
    
    #Adding filter array to DataFrame
    ECO.loc[:, 'filters'] = filters
    
    ### Exposure time check
    exptime_arr = [9399, 3671, 3331, 1319, 2055, 2236, 1758, 10337, 2045, 1237, 
                   2290, 3853, 1928101311, 73024829, 275363, 1241, 31705, 
                   26575,6021, 3548, 3723, 2053, 2249, 3368, 5275, 4069, 
                   171413, 31062, 11431, 5789, 8520, 10071, 6677, 24445, 12605,
                   10757, 50294]
    exp_fil_dict = dict(zip(filters_unique, exptime_arr ))
    
    contents = pd.read_csv(revtxt,header=None,names=['filename'])
    contents.filename = 'http://hla.stsci.edu/cgi-bin/' + contents.filename\
    .astype(str)
    
    #Match and return all columns associated with this ECOID and filename
    #from ECO catalog
    ECO2 = ECO.loc[(ECO.filename.isin(contents.filename)) & \
                   (ECO.ECOID==objname),: ] 
    
    ECOID_groups = ECO2.groupby('filters') 
    ECO_keys     = ECOID_groups.groups.keys() 
     
    ECO_match3 = []
    for key in ECO_keys:
        if ECOID_groups.get_group(key).exptime.sum() >= exp_fil_dict[key]: 
            ECO_match3.append(key) #"good" keys
    ECO_match3 = np.array(ECO_match3)
    
    ### At least 2 filter check
    if len(ECO_match3) >= 2:
        result = True
        with open(path_to_interim + 'goodObj.txt', 'a') as newfile:
            newfile.write(np.unique(ECO2.ECOID)[0]+'\n')
    else:
        result = False
        with open(path_to_interim + 'badObj.txt', 'a') as newfile:
            newfile.write(np.unique(ECO2.ECOID)[0]+'\n')
    
    return result
    
