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

goodObj = '../../../data/interim/goodObj.txt'
percent_blank = '../../../data/interim/percent_blankv2.txt'

#Read goodObj.txt
contents = pd.read_csv(goodObj,header=None)
contents.columns = ['ECO_ID']
#An array of the ECOIDs obtained from goodObj.txt
objs_arr = contents.ECO_ID.values

#Read percent_blank.txt
contents = pd.read_csv(percent_blank, header=None,sep=',')
contents.columns = ['ECO_ID','filter','percent_blank']
contents['percent_blank'] = contents['percent_blank'].map(lambda x: x.rstrip('%'))
contents['percent_blank'] = pd.to_numeric(contents['percent_blank'],errors='coerce')

#Group contents by ECO_ID and filter i.e per coadd per filter
pblank_groups = contents.groupby(['ECO_ID','filter'])
pblank_keys = sorted(pblank_groups.groups.keys())

#Get all keys that correspond to that coadd being 100% blank
pblank100_keys = []
for key_i in pblank_keys:
    pblank = pblank_groups.get_group(key_i)['percent_blank'].values[0]
    if pblank == 100:
        pblank100_keys.append(key_i)

#Iterate through all object folders and multiply the coadds for that object 
#only if it is not 100% blank
for index,obj in enumerate(objs_arr):
    print('Object {0}/110'.format(index+1))
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir(obj)
    print('Getting all coadds')
    coadds_arr = glob('*coadd.fits')
    
    #For a coadd check if it matches with the keys of 100% blank coadds
    print('Checking if coadd is 100% blank')
    counter = 0
    coadds_arr2 = []
    for coadd in coadds_arr:
        blank = False
        for key in pblank100_keys:
            if key[0] in coadd and key[1] in coadd:
                counter+=1
                blank = True
                break
        if not blank:
            coadds_arr2.append(coadd)    
        
        
        
    print(coadds_arr2)
    
    #If after this check only one coadd remains then that means that object has 
    #imaging in less than one filter and it is no longer usable. This should be
    #11.
    if len(coadds_arr2) <=1:
        print('{0} can no longer be used'.format(obj))
        os.chdir('..')
        with open('percent_blank_combinedcoaddv3.txt','a') as newfile:
            newfile.write('{0} can no longer be used\n'.format(obj))
            newfile.close()
           
    #Multiply only if there are 2 or more coadds after this check.
    if len(coadds_arr2) >= 2:
        print('Multiplying coadds')
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu1 = fits.open(coadds_arr2[0])
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        hdu2 = fits.open(coadds_arr2[1])
        sci1 = hdu1[0].data
        sci2 = hdu2[0].data
        
        scidata_temp = sci1*sci2
    
        if len(coadds_arr2) > 2:
            for coadd in coadds_arr2[2:]:
                warnings.simplefilter('ignore', category=AstropyUserWarning)
                hdu = fits.open(coadd)
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
        with open('percent_blank_combinedcoaddv3.txt','a') as newfile:
            newfile.write('{0},{1}\n'.format(obj,percentage))
            newfile.close()

print('Number of coadds that are all blank: {0}'.format(counter))