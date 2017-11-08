#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 12:16:28 2017

@author: asadm2
"""
import pandas as pd
from astropy.io import fits
from glob import glob
import os
import numpy as np
from astropy.wcs import WCS

ECOnew  = '../../../data/interim/ECO_formatted.txt'
goodObj = '../../../data/interim/goodObj.txt'

ECO_new  = pd.read_csv(ECOnew, delimiter='\t')
    
good_Obj         = pd.read_csv(goodObj,header=None)
good_Obj.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
arr_goodObj             = good_Obj.ECO_ID.values

for index,obj in enumerate(arr_goodObj):
    print('Object {0}/{1}'.format(index+1,len(arr_goodObj)))
    print(obj)
    dir_path = os.getcwd()
    #Since this file is in src/ and we need to access the images we have to 
    #move to the data/interim folder the first time this code runs
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    #After the first time all that has to be done is to move out of the  
    #interim ECO object dir and into the next interim ECO object dir
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    
    arr_imgs = glob('*_test_*')
    
    RA = np.unique(ECO_new.RA.loc[ECO_new.ECOID==obj])[0]
    DEC = np.unique(ECO_new.DEC.loc[ECO_new.ECOID==obj])[0]
    
    print('Doing image check')
    good_img_counter = 0
    for filename in arr_imgs:
        len_original = len(arr_imgs)
        hdulist = fits.open(filename)
        header = hdulist[1].header
        w = WCS(header)
        px,py = w.wcs_world2pix(RA,DEC,1)
        px = int(px)
        py = int(py)
        if px <= header['NAXIS1'] and px >= 0 and py <= header['NAXIS2'] and \
        py >= 0:
            good_img_counter += 1
        hdulist.close()
    
    os.chdir('../')
    print('Writing results to text file')
    with open('Obj_in_Img_results.txt', 'a') as newfile:
        newfile.write(obj+' :{0}/{1}\n'.format(good_img_counter,len_original))
    newfile.close()
        