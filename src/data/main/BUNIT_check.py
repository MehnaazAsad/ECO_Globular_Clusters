#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 16:46:35 2017

@author: asadm2
"""

import pandas as pd
import os
from astropy.io import fits
from glob import glob
import warnings
from astropy.utils.exceptions import AstropyUserWarning

goodObj = '../../../data/interim/goodObj.txt'
#Read goodObj.txt
filename1 = pd.read_csv(goodObj,header=None)
filename1.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
objs_arr = filename1.ECO_ID.values

for index,obj in enumerate(objs_arr):
    print('Object {0}/110'.format(index+1))
    print(obj)
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/raw/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir('../'+obj)
    
    arr = glob('*_*')

    for image in arr:
        try:
            warnings.simplefilter('ignore', category=AstropyUserWarning)
            hdu = fits.open(image,ignore_missing_end = True)
            if len(hdu) == 1:
                hdr = hdu[0].header
            else:
                hdr = hdu[1].header
            if ['BUNIT'] in hdr:
                bunit = hdr['BUNIT']
            else:
                bunit = hdr['D001OUUN']
            os.chdir('/fs1/masad/Research/Repositories/ECO_Globular_Clusters/data/interim/')
            with open('BUNITresults.txt', 'a') as filename:
                filename.write(bunit + '\n')
            hdu.close()
            os.chdir('../raw/'+obj)
        except IOError as e:
            pass
        
    filename.close()