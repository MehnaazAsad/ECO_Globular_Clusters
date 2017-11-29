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

goodObj = '../../../data/interim/goodObjv2.txt'
good_Obj = pd.read_csv(goodObj,header=None)
good_Obj.columns = ['ECOID']

ECOnew  = '../../../data/interim/ECO_formatted.txt'
ECO_new  = pd.read_csv(ECOnew, delimiter='\t')

ECOphotcat = '../../../data/interim/ECOphot_final.txt'
ECO_phot_cat  = pd.read_csv(ECOphotcat, delim_whitespace=True,header=None,\
                            skiprows=1,names=['ECOID','RA','DEC','VLG','umag',\
                                              'gmag','rmag','imag','zmag',\
                                              'Jmag','Hmag','Kmag'])

good_Obj_subset = good_Obj[:9]['ECOID'].values

sdssr_calc = []
sdssr_cat = []
y_err = []
for index,obj in enumerate(good_Obj_subset):
    print('{0}/{1} {2}'.format(index+1,len(good_Obj_subset),obj))
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/interim/'+obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir('../'+obj)
        
    comb_coadd = glob(obj+'_comb_coadd.fits')[0]
    f814_coadd = glob(obj+'_acs_wfc_f814w_coadd.fits')[0]
    
    print('Starting source extractor')
    subprocess.call(['sex',comb_coadd+","+f814_coadd,'-ANALYSIS_THRESH','1.5',\
    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1','-PIXEL_SCALE','0.0',\
    '-CATALOG_NAME',obj+'_acs_wfc_f814w.cat']) #CHANGE THIS
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
    print(wcs)
    xx,yy = wcs.wcs_world2pix(ra,dec,1)
    print(xx,yy)

    print('Reading test.cat')
    #CHANGE NAME
    f814w_cat = pd.read_csv(obj+'_acs_wfc_f814w.cat',header=None,\
                            delim_whitespace=True,\
                            names=['petro_mag','petro_magerr','petro_radius',\
                                   'xmin_image','ymin_image','xmax_image',\
                                   'ymax_image','x_image','y_image','x_world',\
                                   'y_world','a_image','class_star'], \
                                   comment='#')
    
#    f814w_cat = pd.read_csv(obj+'_acs_wfc_f814w.cat',header=None,\
#                            delim_whitespace=True,skiprows=7,\
#                            names=['petro_mag','petro_magerr','petro_radius',\
#                                   'xmin_image','ymin_image','xmax_image',\
#                                   'ymax_image','x_image','y_image','a_image',\
#                                   'class_star'])
#    
#    f814w_cat['xmin_image'] = pd.to_numeric(f814w_cat['xmin_image'],errors='coerce')
#    f814w_cat['ymin_image'] = pd.to_numeric(f814w_cat['ymin_image'],errors='coerce')
#    f814w_cat['xmax_image'] = pd.to_numeric(f814w_cat['xmax_image'],errors='coerce')
#    f814w_cat['ymax_image'] = pd.to_numeric(f814w_cat['ymax_image'],errors='coerce')
    
    print('Getting petro mag from test.cat')
#    print(f814w_cat.x_image==xx)
#    print(f814w_cat.y_image==yy)
    f814mag = f814w_cat.petro_mag.loc[((f814w_cat.xmin_image < [xx])&([xx] < \
                                       f814w_cat.xmax_image))&\
                                      ((f814w_cat.ymin_image < [yy])&([yy] < \
                                       f814w_cat.ymax_image))]\
                                      .values[0] 
                                      
    magerr = f814w_cat.petro_magerr.loc[((f814w_cat.xmin_image < [xx])&([xx] <\
                                       f814w_cat.xmax_image))&\
                                      ((f814w_cat.ymin_image < [yy])&([yy] < \
                                       f814w_cat.ymax_image))]\
                                      .values[0]
    f814mag = pd.to_numeric(f814mag)
    print(f814mag)
    print(magerr)
    
    print('Calculating rmag')
    f814mag += zpt814
    sdss_r = f814mag - 1
    
    print('Retrieving rmag from catalog')
    sdss_r_cat = ECO_phot_cat.rmag.loc[ECO_phot_cat.ECOID==obj].values[0]
    
    sdssr_calc.append(sdss_r)
    sdssr_cat.append(sdss_r_cat)
    y_err.append(magerr)
    hdu_f814w_coadd.close()

os.chdir('..')
print('Plotting')
x = np.linspace(0,len(good_Obj_subset)+1,len(good_Obj_subset))
my_xticks = good_Obj_subset
fig1 = plt.figure(figsize=(10,8))
plt.xticks(x, my_xticks,rotation=90)
plt.scatter(x,sdssr_calc, c='r',label='calculated rmag')
plt.scatter(x,sdssr_cat, c='g',label='catalog rmag')
plt.xlabel('ECOID')
plt.ylabel('rmag')
plt.gca().invert_yaxis()
plt.legend()
plt.title('Comparison of calculated rmag and rmag from ECO photometric table')
plt.savefig('calcrmag_catrmag.png')
    
    
    