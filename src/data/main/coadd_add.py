#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:05:26 2017

@author: asadm2
"""

### DESCRIPTION
#This script excludes all coadds that were found to be 100% blank from the
#text file that swarp.py's percent_blank function outputs before adding
#the coadds to make a combined coadd i.e. all the coadds for an object put 
#together. It also exlcudes the list of objects that can no longer be used
#because enough coadds were 100% blank that there was no longer imaging left
#in at least 2 filters. This bad_obj array comes from the text file that 
#coadd_multiply.py output.

from __future__ import division, print_function, absolute_import
from astropy.utils.exceptions import AstropyUserWarning
from astropy.io import fits
from astropy import wcs
from glob import glob
import pandas as pd
import warnings
import os

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'

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

#Object IDs from percent_blank_combined_coaddv3.txt for objects that have 
#enough coadds that are 100% blank that the unique filter count drops below 2
#(Probably better to have an automated way to extract this list?)
bad_obj = ['ECO03777', 'ECO03847', 'ECO04454', 'ECO06176', 'ECO06184', \
           'ECO06685', 'ECO07159', 'ECO08944', 'ECO09176', 'ECO10556', \
           'ECO12253']

#Iterate through all object folders and add the coadds for that object 
#only if it is not 100% blank
for index,obj in enumerate(arr_goodObj): 
    print('Object {0}/{1}'.format(index+1,len(arr_goodObj)))
    if obj in bad_obj:
        print('Bad object: {0}'.format(obj))
        continue
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir(path_to_interim + obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    print('Getting all coadds')
    
    #Remove older version of this coadd before it is made 
    os.remove(obj+'comb_coadd.fits')
    
    #Grab all coadds in current folder
    coadds_arr = glob('*coadd.fits')
    
    #For a coadd check if it matches with the keys of 100% blank coadds. If it
    #does then increment the counter otherwise apped the coadd to arr_2.
    print('Checking if coadd is 100% blank')
    coadds_arr2 = []
    for coadd in coadds_arr:
        blank = False
        for key in pblank100_keys:
            if key[0] in coadd and key[1] in coadd:
                blank = True
                break
        if not blank:
            coadds_arr2.append(coadd) 
            
    #Combine images (add) only if there are 2 or more coadds after this check.
    if len(coadds_arr2) >= 2:      
        print('Adding coadds')
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu1 = fits.open(coadds_arr2[0])
        w = wcs.WCS(hdu1[0].header)
        header = w.to_header()
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu2 = fits.open(coadds_arr2[1])
        sci1 = hdu1[0].data
        sci2 = hdu2[0].data
        
        #Add the first two coadds first
        scidata_temp = sci1+sci2
        
        #If there were more than two then add scidata_temp with each
        #consecutive image
        if len(coadds_arr2) > 2:
            for coadd in coadds_arr2[2:]:
                warnings.simplefilter('ignore', category=AstropyUserWarning)
                hdu = fits.open(coadd)
                scidata = hdu[0].data
                scidata_temp += scidata
                hdu.close()
        
        #Data of combined coadd is the average of all the individual coadds 
        scidata_temp = scidata_temp/(len(coadds_arr2))
        
        #Write to new fits image
        fits.writeto(obj+'_comb_coadd.fits', scidata_temp, header,\
                     output_verify='ignore')
                
        hdu1.close()
        hdu2.close()
            
    os.chdir('..')