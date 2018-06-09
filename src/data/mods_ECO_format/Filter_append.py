# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:53:32 2017

@author: asadm2
"""

import pandas as pd
import numpy as np

def filter_append(filename):
    """
    This function takes in the original ECO catalog and reformats it to 
    include a column for the filter that was used to take the image.
    
    Args: 
        filename: Original ECO catalog
    
    Returns:
        Newly formatted ECO catalog
    """
    global ECO
    ECO = pd.read_csv(filename, delimiter='\s+', header=None, \
                      names=['ECOID', 'HSTOBJ', 'RA', 'DEC', 'exptime', \
                             'camera', 'filename'])
    
    ECO['exptime'] = pd.to_numeric(ECO['exptime'],errors='coerce')
    ECO['filename'] = ECO['filename'].astype('str')
    
    # Number of elements in ECO
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
    
    #Adding filter array to DataFrame
    ECO.loc[:, 'filters'] = filters
    
    return ECO