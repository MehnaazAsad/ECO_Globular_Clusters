# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 16:39:59 2017

@author: asadm2
"""
from __future__ import division, print_function, absolute_import
import astropy
from glob import glob
from astropy.io import fits
from astropy.io.fits import getdata
import warnings
from astropy.utils.exceptions import AstropyUserWarning
import pandas as pd
import numpy as np
import os
import subprocess
import imp

foo = imp.load_source('ECO_funcs', '../mods_swarp/ECO_funcs.py')

#Formatted ECO catalog with filters and new_filename column and list of 
#good ECO objects
ECOnew  = '../../../data/interim/ECO_formatted.txt'
goodObj = '../../../data/interim/goodObjv2.txt'

#Using the make_dict function to return a dictionary whose keys are ECOID and
#filter name and whose values are filenames, return the list of sorted keys,
#and the formatted ECO catalog
ECO_dict,ECO_keys,ECOcat = foo.make_dict(ECOnew,goodObj)

for index,key in enumerate(ECO_keys):
    print('Key {0}/{1}'.format(index+1,len(ECO_keys)))
    print(key)
    obj = key[0]
    filter_i = key[1]
    dir_path = os.getcwd()
    #Since this file is in src/ and we need to access the images we have to 
    #move to the data/raw folder the first time this code runs
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/raw/'+obj)
    #After the first time all that has to be done is to move out of the  
    #interim ECO object dir and into the raw ECO object dir
    elif os.path.basename(dir_path) != obj:
        print(dir_path)
        os.chdir('../raw/'+obj)
    #Using the ECOID,filter key to acess the values in the new_filename column
    #i.e. get a list of the images associated with this key pair
    print("Getting images from catalog")
    imgs_from_cat = ECO_dict.get_group(key)['new_filename'].values

    #THIS WON'T WORK AS THEY HAVE TO BE SEPARATED BY FILTER. KEEP LINE 49 and
    #next 3 lines and have a loop afterwards that finds common filenames in 
    #both arrays and appends to new array that will be used as the imgs array
    #that update_header takes. If no match then continue to next loop iteration
    print("Getting images from text file")
    goodimagesfile = obj+'_goodimages.txt'
    goodimages = pd.read_csv(goodimagesfile,header=None,names=['filenames'])
    imgs_from_txt = goodimages.values
    
    print("Common images are")
    imgs = []
    for img in imgs_from_cat:
        if img in imgs_from_txt:
            print(img)
            img = img + '.fits'
            imgs.append(img)
        
    #Using the update_header function to output the name of the text file
    #that contains the names of the new versions of the fits images being
    #input into the function
    img_list_txt = foo.update_header(imgs,obj,filter_i)

    #To check for case where all images were in Except clause of update_header
    #and run swarp even for single images
    if len(imgs) >= 1 and img_list_txt != 'error': 
        ra = np.unique(ECOcat.RA.loc[ECOcat.ECOID==obj])[0]
        dec = np.unique(ECOcat.DEC.loc[ECOcat.ECOID==obj])[0]
        
        #Running swarp with parameters. The only ones that change are the list
        #of images, the RA and DEC of the object, the object name and the 
        #filter name
        print('Starting to swarp')    
        subprocess.call(['swarp','@'+img_list_txt,'-CENTER_TYPE','MANUAL', \
        '-CENTER',str(ra)+","+str(dec),'-NTHREADS','32', '-COPY_KEYWORDS', \
        'PHOTFLAM,PHOTZPT,PHOTPLAM','-FSCALASTRO_TYPE','VARIABLE',\
        '-COMBINE_BUFSIZE','10000','-MEM_MAX','1000','-VMEM_MAX','5000',\
        '-IMAGE_SIZE','3960','-PIXELSCALE_TYPE','MANUAL','-PIXEL_SCALE',\
        '0.05','-GAIN_KEYWORD','CCDGAIN','-GAIN_DEFAULT','1.0',\
        '-SUBTRACT_BACK','Y','-IMAGEOUT_NAME',\
        obj+'_'+filter_i+'_'+'coadd.fits','-WEIGHTOUT_NAME',\
        obj+'_'+filter_i+'_'+'coadd.weight.fits','-XML_NAME',\
        obj+'_'+filter_i+'_'+'swarp.xml'])
        print('Finished swarping')
        
        coadd = fits.open(obj+'_'+filter_i+'_'+'coadd.fits')
        #Open the coadd and use the percent_blank function 
        foo.percent_blank(coadd,obj,filter_i)
        #So that you can distinguish between interim object folder and raw 
        #object filder
        
    else:
        os.chdir('..') # Moving one level up when no imgs are found
        
    #If the length of imgs is 1 then check if the new version that 
    #update_header function creates exists. If it does then calculate blank
    #percentage. If the length is 0 i.e. since the original ECO catalog is 
    #being used some ECOID, filter pairs no longer exist since they were
    #eliminated in preliminary exptime and >2filter checks so for those pairs
    #the files won't exist. They are ignored which is why percent_blank doesn't
    #print anything to the terminal.
#    else:
#        filename = imgs[0]+'_test_'+filter_i+'.fits'
#        if os.path.exists(filename): #avoid empty cases where files don't exist 
#            single_img = fits.open(filename)
#            foo.percent_blank(single_img,obj,filter_i)       