#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 23:25:33 2017

@author: asadm2
"""

from __future__ import division, print_function, absolute_import
from glob import glob
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import subprocess

eco_repo_path = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'

goodObj = '../../../data/interim/goodObj.txt'
good_Obj = pd.read_csv(goodObj,header=None)
good_Obj.columns = ['ECOID']

ECOnew  = '../../../data/interim/ECO_formatted.txt'
ECO_new  = pd.read_csv(ECOnew, delimiter='\t')

ECOphotcat = '../../../data/interim/ECOphot_final.txt'
ECO_phot_cat  = pd.read_csv(ECOphotcat, delim_whitespace=True,header=None,\
                            skiprows=1,names=['ECOID','RA','DEC','VLG','umag',\
                                              'gmag','rmag','imag','zmag',\
                                              'Jmag','Hmag','Kmag'])

good_Obj_new = good_Obj.ECOID[:4].append(good_Obj.ECOID[5:11].append(good_Obj\
                             .ECOID[104:105]))

sdssr_calc = []
sdssr_cat = []
for obj in good_Obj_new:
    print(obj)
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir('../interim/'+obj)
        
    comb_coadd = glob(obj+'comb_coadd.fits')[0]
    f814_coadd = glob(obj+'_acs_wfc_f814w_coadd.fits')[0]
    
    print('Starting source extractor')
    subprocess.call(['sex',comb_coadd+","+f814_coadd,'-ANALYSIS_THRESH','1.5',\
    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1',\
    '-CATALOG_NAME',obj+'acs_wfc_f814w.cat'])
    print('Finished running source extractor')
    
    hdu_f814w_coadd = fits.open(f814_coadd)  
    prihdr = hdu_f814w_coadd[0].header
    photflam814 = prihdr['PHOTFLAM']
    photzpt814 = prihdr['PHOTZPT']
    print('Calculating zeropoint')
    zpt814 = -2.5*np.log10(photflam814) + photzpt814

    ra = np.unique(ECO_new.RA.loc[ECO_new.ECOID==obj])[0]
    dec = np.unique(ECO_new.DEC.loc[ECO_new.ECOID==obj])[0]
    
    print('Converting to pixel coordinates')
    wcs = WCS(prihdr)   
    xx,yy = wcs.all_world2pix(ra,dec,1)
    
#remember to copy already uncommented default.param file and remove other one
#remember to copy all default files to last obj in contents_new

#glob comb_coadd and store as image1
#glob f814 coadd and store as image2
#use subprocess to pass parameters
#open filter image and get necessary keywords for finding zeropoint
#get RA and DEC of galaxy from ECO catalog
#open filter image and get WCS info to find world2pix transformation for RA and
#DEC of galaxy
#open catalog that sextractor generated and get petro mag in f814 at 
#calculated x and y pixel values
#Apply zeropoint to magnitude
#Use conversion to change to r band
#Get r band mag from photometry catalog
#Plot both points for each galaxy

#need to use WCS

    print('Reading test.cat')
    f814w_cat = pd.read_csv(obj+'acs_wfc_f814w.cat',header=None,\
                            delim_whitespace=True,skiprows=1,\
                            names=['petro_mag','petro_magerr','petro_radius',\
                                  'x_image','y_image','a_image','class_star'])
    
    print('Getting petro mag from test.cat')
    f814mag = f814w_cat.petro_mag.loc[(f814w_cat.x_image==xx)&\
                                      (f814w_cat.y_image==yy)].values[0]
    
    print('Calculating rmag')
    f814mag += zpt814
    sdss_r = f814mag - 1
    
    print('Retrieving rmag from catalog')
    sdss_r_cat = ECO_phot_cat.rmag.loc[ECO_phot_cat.ECOID==obj].values[0]
    
    sdssr_calc.append(sdss_r)
    sdssr_cat.append(sdss_r_cat)

print('Plotting')
x = np.linspace(0,len(good_Obj_new.values)+1,len(good_Obj_new.values))
my_xticks = good_Obj_new.values
fig1 = plt.figure(figsize=(10,8))
plt.xticks(x, my_xticks,rotation=90)
plt.scatter(x,sdssr_calc, c='r',label='calculated rmag')
plt.scatter(x,sdssr_cat, c='g',label='catalog rmag')
plt.xlabel('ECOID')
plt.ylabel('rmag')
plt.gca().invert_yaxis()
plt.title('Comparison of calculated rmag and rmag from ECO photometric table')
plt.savefig('calcrmag_catrmag.png')
    
    
    