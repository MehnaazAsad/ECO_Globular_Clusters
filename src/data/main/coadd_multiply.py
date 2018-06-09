#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 00:05:59 2017

@author: asadm2
"""
### DESCRIPTION
#This script excludes all coadds that were found to be 100% blank from the
#text file that swarp.py's percent_blank function outputs before multiplying
#the coadds to calculate a new percentage for the combined coadd i.e. all the
#coadds for an object put together. This result is output as a text file as 
#well. However, if after accounting for those coadds that were 100% blank 
#there aren't coadds in atleast 2 filters then that object is completely 
#discarded and should be removed from the goodObjv2 text file that has been
#used until this point.

from __future__ import division, print_function, absolute_import
from astropy.utils.exceptions import AstropyUserWarning
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

#If re-running script then this file has to be removed since
#it is being opened in 'append' mode in this script
if os.path.isfile(path_to_interim + 'percent_blank_combinedcoaddv3.txt'): 
    os.remove(path_to_interim + 'percent_blank_combinedcoaddv3.txt')

goodObj = path_to_interim + 'goodObjv2.txt'
percentblank = path_to_interim + 'percent_blankv2.txt'

#Read goodObj.txt
good_Obj = pd.read_csv(goodObj,header=None, names = ['ECO_ID'])
#An array of the ECOIDs obtained from goodObj.txt
arr_goodObj = good_Obj.ECO_ID.values

#Read percent_blankv2.txt
percent_blank = pd.read_csv(percentblank, header=None,\
                            names=['ECO_ID','filter','percent_blank'],sep=',')
percent_blank['percent_blank'] = percent_blank['percent_blank'].map(lambda x: \
             x.rstrip('%'))
percent_blank['percent_blank'] = pd.to_numeric(percent_blank['percent_blank'],\
             errors='coerce')

#Group contents by ECO_ID and filter i.e per coadd per filter
pblank_groups = percent_blank.groupby(['ECO_ID','filter'])
pblank_keys = sorted(pblank_groups.groups.keys())

#Get all keys that correspond to that coadd being 100% blank
pblank100_keys = []
for key_i in pblank_keys:
    pblank = pblank_groups.get_group(key_i)['percent_blank'].values[0]
    if pblank == 100:
        pblank100_keys.append(key_i)

#Iterate through all object folders and multiply the coadds for that object 
#only if it is not 100% blank
for index,obj in enumerate(arr_goodObj):
    print('Object {0}/{1}'.format(index+1,len(arr_goodObj)))
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir(path_to_interim + obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    print('Getting all coadds')
    coadds_arr = glob('*coadd.fits')
    
    #For a coadd, check if it matches with the keys of 100% blank coadds. If it
    #does then increment the counter otherwise apped the coadd to arr_2 to be
    #used later.
    print('Checking if coadd is 100% blank')
    counter = 0
    coadds_arr_2 = []
    for coadd in coadds_arr:
        blank = False
        for key in pblank100_keys:
            if key[0] in coadd and key[1] in coadd:
                counter+=1
                blank = True
                break
        if not blank:
            coadds_arr_2.append(coadd)    
        
    #If after this check only one coadd remains then that means that object has 
    #imaging in less than one filter and it is no longer usable. 
    if len(coadds_arr_2) <=1:
        print('{0} can no longer be used'.format(obj))
        os.chdir(path_to_interim)
        with open('percent_blank_combinedcoaddv3.txt','a') as newfile:
            newfile.write('{0} can no longer be used\n'.format(obj))
            newfile.close()
           
    #Multiply only if there are 2 or more coadds after this check.
    if len(coadds_arr_2) >= 2:
        print('Multiplying coadds')
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu1 = fits.open(coadds_arr_2[0])
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu2 = fits.open(coadds_arr_2[1])
        sci1 = hdu1[0].data
        sci2 = hdu2[0].data
        
        #Multiply the first two coadds first
        scidata_temp = sci1*sci2
        
        #If there were more than two then multiply scidata_temp with each
        #consecutive image
        if len(coadds_arr_2) > 2:
            #sliced since first 2 images have already been multiplied above
            for coadd in coadds_arr_2[2:]: 
                warnings.simplefilter('ignore', category=AstropyUserWarning)
                hdu = fits.open(coadd)
                scidata = hdu[0].data
                scidata_temp *= scidata
                hdu.close()
            
        hdu1.close()
        hdu2.close()
        
        #Calculating percentage blank in the combined image of all coadds
        new  = np.ravel(scidata_temp) #returns 1-D array
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
            #Go here if there aren't any pixel values that are exactly 0 in 
            #which case the image is 0% blank
            percentage = 0
        
        print('Percentage blank {0}%'.format(percentage))
        
        os.chdir(path_to_interim)      
        with open('percent_blank_combinedcoaddv3.txt','a') as newfile:
            newfile.write('{0},{1}\n'.format(obj,percentage))
            newfile.close()

    print('Number of coadds that are all blank: {0}'.format(counter))
