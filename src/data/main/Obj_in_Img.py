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
import warnings
from astropy.utils.exceptions import AstropyUserWarning

ECOnew  = '../../../data/interim/ECO_formatted.txt'
goodObj = '../../../data/interim/goodObj.txt'

ECO_new  = pd.read_csv(ECOnew, delimiter='\t')  
    
good_Obj         = pd.read_csv(goodObj,header=None)
good_Obj.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
arr_goodObj             = good_Obj.ECO_ID.values

exptime_arr = [9399, 3671, 3331, 1319, 2055, 2236, 1758, 10337, 2045, 1237, 
                   2290, 3853, 1928101311, 73024829, 275363, 1241, 31705, 26575,
                   6021, 3548, 3723, 2053, 2249, 3368, 5275, 4069, 171413, 31062,
                   11431, 5789, 8520, 10071, 6677, 24445, 12605, 10757, 50294]
filters_unique = np.unique(ECO_new.filters.values)

exp_fil_dict = dict(zip(filters_unique, exptime_arr ))

arr_good_img_counter = []
for index,obj in enumerate(arr_goodObj):
    print('Object {0}/{1}'.format(index+1,len(arr_goodObj)))
    print(obj)
    dir_path = os.getcwd()
    #Since this file is in src/ and we need to access the images we have to 
    #move to the data/raw folder the first time this code runs
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/raw/'+obj)
    #After the first time all that has to be done is to move out of the  
    #interim ECO object dir and into the next interim ECO object dir
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    
    if os.path.isfile(obj+'_goodimages.txt'): 
        os.remove(obj+'_goodimages.txt')
    
    arr_imgs = glob('*.fits')
    
    RA = np.unique(ECO_new.RA.loc[ECO_new.ECOID==obj])[0]
    DEC = np.unique(ECO_new.DEC.loc[ECO_new.ECOID==obj])[0]
    
    print('Doing image check')
    good_img_counter = 0
    good_img_arr = []
    for image in arr_imgs:
        len_original = len(arr_imgs)
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        try:
            hdulist = fits.open(image,ignore_missing_end=True)
            header = hdulist[1].header
            w = WCS(header)
            px,py = w.wcs_world2pix(RA,DEC,1)
            px = int(px)
            py = int(py)
            if px <= header['NAXIS1'] and px >= 0 and py <= header['NAXIS2'] and \
            py >= 0:
                scidata = hdulist[1].data
                val_at_pix = scidata[py-1,px-1]
                if val_at_pix != 0:
                    good_img_counter += 1
                    os.chdir('../../interim/'+obj)
                    new_filename = os.path.splitext(image)[0]
                    with open(obj+'_goodimages.txt','a') as newfile:
                        newfile.write(new_filename+'\n')
                    newfile.close()
                    os.chdir('../../raw/'+obj)
                    good_img_arr.append(new_filename)
            hdulist.close()
        except IOError as e:
            dir_path = os.getcwd()
            if os.path.basename(dir_path) == obj:
                os.chdir('../../interim/')
            with open('Error_objinimg.txt','a') as newfile:              
                newfile.write('Object {0} and image {1} raises {2} error.\n'.format\
                (obj,image,e))
            newfile.close()
            os.chdir('../raw/'+obj)
    
    arr_good_img_counter.append(good_img_counter)    
    
    os.chdir('../../interim')
    print('Writing results to text file')
    with open('Obj_in_Img_results.txt', 'a') as newfile:
        newfile.write(obj+' {0} {1}\n'.format(good_img_counter,len_original))
    newfile.close()
    
    ECO2 = ECO_new.loc[(ECO_new.new_filename.isin(good_img_arr)) & \
                       (ECO_new.ECOID==obj),: ]
    ECOID_groups = ECO2.groupby('filters') #returns dict of grouped dataframes
    ECO_keys     = ECOID_groups.groups.keys() #get the key list of each group in the dict
    
    ECO_match = []
    final_good_img_arr = []
    for key in ECO_keys:
        if ECOID_groups.get_group(key).exptime.sum() >= exp_fil_dict[key]: #exptime check
            ECO_match.append(key) #"good" keys
            final_good_img_arr.append(ECOID_groups.get_group(key).new_filename)
        
    final_good_img_num = sum(len(x) for x in final_good_img_arr)
    
    if len(ECO_match) >= 2:
        filter_num = len(ECO_match)
        with open('Expfil2_results_good.txt', 'a') as newfile:
            newfile.write(obj+' {0} {1} {2}\n'.format(final_good_img_num,\
                          len_original,filter_num))
    else:
        with open('Expfil2_results_bad.txt', 'a') as newfile:
            newfile.write(obj+'\n')
    
    os.chdir('../raw/')
    
    