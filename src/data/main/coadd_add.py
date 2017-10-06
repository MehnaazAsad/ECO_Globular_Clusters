#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 21:05:26 2017

@author: asadm2
"""

from __future__ import division, print_function, absolute_import
import pandas as pd
import os
from glob import glob
from astropy.io import fits
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from astropy import wcs

goodObj = '../../../data/interim/goodObj.txt'
percent_blank = '../../../data/interim/percent_blankv2.txt'

#Read goodObj.txt
filename1 = pd.read_csv(goodObj,header=None)
filename1.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
objs_arr = filename1.ECO_ID.values

#Read percent_blank.txt
filename2 = pd.read_csv(percent_blank, header=None,sep=',')
filename2.columns = ['ECO_ID','filter','percent_blank']
filename2['percent_blank'] = filename2['percent_blank'].map(lambda x: x.rstrip('%'))
filename2['percent_blank'] = pd.to_numeric(filename2['percent_blank'],errors='coerce')

#Group contents by ECO_ID and filter i.e per coadd per filter
pblank_groups = filename2.groupby(['ECO_ID','filter'])
pblank_keys = sorted(pblank_groups.groups.keys())

#Get all keys that correspond to that coadd being 100% blank
pblank100_keys = []
for key_i in pblank_keys:
    pblank = pblank_groups.get_group(key_i)['percent_blank'].values[0]
    if pblank == 100:
        pblank100_keys.append(key_i)

#Object IDs from percent_blank_combined_coaddv2.txt for objects that have 
#enough coadds that are 100% blank that the unique filter count drops below 2
bad_obj = ['ECO03777', 'ECO03847', 'ECO04454', 'ECO06176', 'ECO06184', \
           'ECO06685', 'ECO07159', 'ECO08944', 'ECO09176', 'ECO10556', \
           'ECO12253']

#Iterate through all object folders and add the coadds for that object 
#only if it is not 100% blank
for index,obj in enumerate(objs_arr[0:1]):
    print('Object {0}/110'.format(index+1))
    if obj in bad_obj:
        print('Bad object {0}'.format(obj))
        continue
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    print('Getting all coadds')
    
#    os.remove(obj+'comb_coadd.fits') #remove after run is complete
    
    coadds_arr = glob('*coadd.fits')
    
    #For a coadd check if it matches with the keys of 100% blank coadds
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
        
        scidata_temp = sci1+sci2
    
        if len(coadds_arr2) > 2:
            for coadd in coadds_arr2[2:]:
                warnings.simplefilter('ignore', category=AstropyUserWarning)
                hdu = fits.open(coadd)
                scidata = hdu[0].data
                scidata_temp += scidata
                hdu.close()
        
        scidata_temp = scidata_temp/(len(coadds_arr2))
            
        fits.writeto(obj+'comb_coadd.fits', scidata_temp, header, output_verify='ignore')
                
        hdu1.close()
        hdu2.close()
            
    os.chdir('..')
        
        
        
        
        
        
        