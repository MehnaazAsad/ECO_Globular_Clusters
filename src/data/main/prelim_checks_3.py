#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 12:16:28 2017

@author: asadm2
"""
### DESCRIPTION
#This script uses the objects that passed prelim_checks_2 (goodObj.txt).
#It repeats the checks carried out in prelim_checks_2 with the additional 
#check of whether the pixel value at the RA and DEC position is 0 or not. 
#This was because for certain images the object was "in" the image but in
#a black region due to the actual image being rotated a certain degree. (Pixel
#value in black region should be 0.) Once these images were filtered out,
#the exposure time check and the "at least 2 filter" check were repeated. The
#remaining "good" objects that passed these checks are written to a text file 
#and the "bad" objects to a separate text file. 

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import WCS
from astropy.io import fits
from glob import glob
import pandas as pd
import numpy as np
import warnings
import os

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'

#If re-running script then these files have to be removed since
#they are being opened in 'append' mode below
if os.path.isfile(path_to_interim + 'Error_prelim3.txt'): 
    os.remove(path_to_interim + 'Error_prelim3.txt')
if os.path.isfile(path_to_interim + 'Prelim3_results.txt'): 
    os.remove(path_to_interim + 'Prelim3_results.txt')
if os.path.isfile(path_to_interim + 'Prelim3_results_good.txt'): 
    os.remove(path_to_interim + 'Prelim3_results_good.txt')
if os.path.isfile(path_to_interim + 'Prelim3_results_bad.txt'): 
    os.remove(path_to_interim + 'Prelim3_results_bad.txt')
if os.path.isfile(path_to_interim + 'goodObjv2.txt'): 
    os.remove(path_to_interim + 'goodObjv2.txt')

#Formatted ECO catalog
ECOnew  = path_to_interim + 'ECO_formatted.txt'
#Text file of objects that passed prelim_checks_2
goodObj = path_to_interim + 'goodObj.txt'
#Reading catalog
ECO_new  = pd.read_csv(ECOnew, delimiter='\t')  
#Reading goodObj.txt   
good_Obj         = pd.read_csv(goodObj,header=None,names=['ECO_ID'])
arr_goodObj      = good_Obj.ECO_ID.values

###Exposure time and filter dictionary
exptime_arr = [9399, 3671, 3331, 1319, 2055, 2236, 1758, 10337, 2045, 1237,               
               2290, 3853, 1928101311, 73024829, 275363, 1241, 31705, 
               26575, 6021, 3548, 3723, 2053, 2249, 3368, 5275, 4069, 
               171413, 31062, 11431, 5789, 8520, 10071, 6677, 24445, 12605,
               10757, 50294]

filters_unique = np.unique(ECO_new.filters.values)
exp_fil_dict = dict(zip(filters_unique, exptime_arr ))
###

arr_good_img_counter = []
for index,obj in enumerate(arr_goodObj):
    
    print('Object {0}/{1}'.format(index+1,len(arr_goodObj)))
    print(obj)
    
    dir_path = os.getcwd()
    
    #Since this script is in src/ and we need to access the images we have to 
    #move to the data/raw folder the first time this code runs
    if os.path.basename(dir_path) == 'main':
        os.chdir(path_to_raw + obj)
    #After the first time all that has to be done is to move out of the  
    #interim ECO object dir and into the next interim ECO object dir
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    
    if os.path.isfile(path_to_interim + obj + '/' + obj + '_goodimages.txt'): 
        os.remove(path_to_interim + obj + '/' + obj + '_goodimages.txt')
    
    #Grab all images in current object directory
    arr_imgs = glob('*.fits')
    
    #Get RA and DEC from catalog
    RA = np.unique(ECO_new.RA.loc[ECO_new.ECOID==obj])[0]
    DEC = np.unique(ECO_new.DEC.loc[ECO_new.ECOID==obj])[0]
    
    ###Object-in-image check
    print('Doing object-in-image check')
    
    good_img_counter = 0
    good_img_arr = []
    for image in arr_imgs:
        len_original = len(arr_imgs)
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        try:
            hdulist = fits.open(image,ignore_missing_end=True)
            #Some images have just one header while others have multiple
            if len(hdulist) == 1:
                header = hdulist[0].header
            else:
                header = hdulist[1].header
            w = WCS(header)
            px,py = w.wcs_world2pix(RA,DEC,1)
            px = int(px)
            py = int(py)
            if px <= header['NAXIS1'] and px >= 0 and py <= header['NAXIS2']\
            and py >= 0:
                if len(hdulist) == 1:
                    scidata = hdulist[0].data
                else:
                    scidata = hdulist[1].data
                #Checking if pixel value is 0 or not
                val_at_pix = scidata[py-1,px-1]
                if val_at_pix != 0:
                    good_img_counter += 1
                    #remove .fits extension
                    new_filename = os.path.splitext(image)[0] 
                    #Write good images to new file that exists in data/raw/
                    #per object
                    with open(obj+'_goodimages.txt','a') as newfile:
                        newfile.write(new_filename+'\n')
                    newfile.close()
                    
                    good_img_arr.append(new_filename)
                    
            hdulist.close()
        #Deal with errors with FITS files
        except IOError as e:
            dir_path = os.getcwd()
            if os.path.basename(dir_path) == obj:
                os.chdir(path_to_interim)
            with open('Error_prelim3.txt','a') as newfile:              
                newfile.write('Object {0} and image {1} raises {2} error.\n'.\
                              format(obj,image,e))
            newfile.close()
            os.chdir(path_to_raw + obj)
    
    arr_good_img_counter.append(good_img_counter)    
    
    os.chdir(path_to_interim)
    
    #Writing ECOID, number of good images AFTER check and total number of 
    #images BEFORE this check to text file.
    print('Writing results to text file')
    with open('Prelim3_results.txt', 'a') as newfile:
        newfile.write(obj+' {0} {1}\n'.format(good_img_counter,len_original))
    newfile.close()
    ###
    
    #Get all columns that match filename and ECOID from ECO catalog
    ECO2 = ECO_new.loc[(ECO_new.new_filename.isin(good_img_arr)) & \
                       (ECO_new.ECOID==obj),: ]
    ECOID_groups = ECO2.groupby('filters') #returns dict of grouped dataframes
    ECO_keys     = ECOID_groups.groups.keys() #get the key list of each group 
                                              #in the dict
    
    ### Exposure time check
    ECO_match = []
    final_good_img_arr = []
    for key in ECO_keys:
        if ECOID_groups.get_group(key).exptime.sum() >= exp_fil_dict[key]: 
            ECO_match.append(key) #"good" keys
            final_good_img_arr.append(ECOID_groups.get_group(key).new_filename)
    ###
    
    #Final number of good images AFTER exposure time check
    final_good_img_num = sum(len(x) for x in final_good_img_arr)
    
    ### "At least 2 filter" check
    if len(ECO_match) >= 2:
        filter_num = len(ECO_match)
        #Writing ECOID, final number of good images AFTER exposure time check,
        #total number of images BEFORE any checks and number of unique filters
        #per object to text file
        with open('Prelim3_results_good.txt', 'a') as newfile:
            newfile.write(obj+' {0} {1} {2}\n'.format(final_good_img_num,\
                          len_original,filter_num))
        newfile.close()
        
        #Writing ECOID of objects that passed all checks to text file
        with open('goodObjv2.txt', 'a') as newfile:
            newfile.write(obj+'\n')
        newfile.close()
            
    else:
        #Writing ECOID of objects that did not pass these checks to text file
        with open('Prelim3_results_bad.txt', 'a') as newfile:
            newfile.write(obj+'\n')
        newfile.close()
    ###
    os.chdir(path_to_raw)
    
    