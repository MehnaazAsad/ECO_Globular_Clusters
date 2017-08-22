#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 00:05:59 2017

@author: asadm2
"""

from __future__ import division, print_function, absolute_import
from glob import glob
from astropy.io import fits
import warnings
from astropy.utils.exceptions import AstropyUserWarning
import pandas as pd
import numpy as np
import os
import imp

foo = imp.load_source('ECO_funcs', '../mods_swarp/ECO_funcs.py')

goodObj = '../../../data/interim/goodObj.txt'

#Read goodObj.txt
ECO_obj = pd.read_csv(goodObj, delimiter='\t') 
contents = pd.read_csv(goodObj,header=None)
contents.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
objs_arr = contents.ECO_ID.values

for obj in objs_arr:
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    coadds_arr = glob('*coadd.fits')
    
    hdu1 = fits.open(coadds_arr[0])
    hdu2 = fits.open(coadds_arr[1])
    sci1 = hdu1[0].data
    sci2 = hdu2[0].data
    
    scidata_temp = sci1*sci2
    
    if len(coadds_arr) > 2:
        for obj in coadds_arr[2:]:
            hdu = fits.open(obj)
            scidata = hdu[0].data
            scidata_temp *= scidata
            hdu.close()
            
    hdu1.close()
    hdu2.close()
    
    new  = np.ravel(scidata_temp)
    new = np.sort(new)
    
    try:
        #Find the range of indices where the pixel value is 0 since the
        #array is sorted and sometimes there are negative values which 
        #we want to include
        index_0_ends = np.where(new == 0)[0][-1]
        index_0_starts = np.where(new == 0)[0][0]
        index_diff = index_0_ends - index_0_starts
        size = len(new)
        #Calculate percentage blank
        percentage = ((index_diff+1)/size)*100
    except:
        #Go here if there aren't any pixel values that are exactly 0 in which 
        #case the image is 0% blank
        percentage = 0
    
    print('Percentage blank {0}%'.format(percentage))
    
    os.chdir('..')
    
    with open('percent_blank_combinedcoadd.txt','a') as newfile:
        newfile.write('{0},{1}\n'.format(obj,percentage))
        newfile.close()
    
    



        

        
        
        
        
    
        
