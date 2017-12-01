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
from astropy import units as u
from astropy.coordinates import SkyCoord

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

good_Obj_subset = good_Obj[:4]['ECOID'].append(good_Obj[5:9]['ECOID']).values

sdssr_petro_calc = []
sdssr_iso_calc = []
sdssr_isocorr_calc = []
sdssr_auto_calc = []
sdssr_cat = []

#y_err = []
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
                            names=['iso_mag','isocorr_mag','auto_mag',\
                                   'petro_flux','petro_fluxerr','petro_mag',\
                                   'petro_magerr','petro_radius',\
                                   'xmin_image','ymin_image','xmax_image',\
                                   'ymax_image','x_image','y_image','x_world',\
                                   'y_world','a_image','class_star'], \
                                   comment='#')
    
    print('Getting petro mag from test.cat')


    cat_sextractor = SkyCoord(f814w_cat['x_world']*u.deg,\
                              f814w_cat['y_world']*u.deg)
    cat_eco = SkyCoord(ra*u.deg, dec*u.deg)
    idx_sdss, d2d_sdss, d3d_sdss = cat_eco.match_to_catalog_sky(cat_sextractor)

    f814mag = f814w_cat.petro_mag.values[idx_sdss]
    isomag = f814w_cat.iso_mag.values[idx_sdss]
    isocorrmag = f814w_cat.isocorr_mag.values[idx_sdss]
    automag = f814w_cat.auto_mag.values[idx_sdss]
    
    petroflux = f814w_cat.petro_flux.values[idx_sdss]
    
    os.chdir('..')
    with open('sextractor_magflux.txt','a') as newfile:
        newfile.write('{0},{1}\n'.format(f814mag,petroflux))
        newfile.close()
    os.chdir(obj)
    
#    f814mag = f814w_cat.petro_mag.loc[((f814w_cat.xmin_image < [xx])&([xx] < \
#                                       f814w_cat.xmax_image))&\
#                                      ((f814w_cat.ymin_image < [yy])&([yy] < \
#                                       f814w_cat.ymax_image))]\
#                                      .values[0] 
#                                      
#    magerr = f814w_cat.petro_magerr.loc[((f814w_cat.xmin_image < [xx])&([xx] <\
#                                       f814w_cat.xmax_image))&\
#                                      ((f814w_cat.ymin_image < [yy])&([yy] < \
#                                       f814w_cat.ymax_image))]\
#                                      .values[0]
#
#    f814mag = pd.to_numeric(f814mag)
    print(f814mag)

    print('Calculating rmag')
    f814mag += zpt814
    sdss_r_petro = f814mag - 1
    
    isomag += zpt814
    sdss_r_iso = isomag - 1
    
    isocorrmag += zpt814
    sdss_r_isocorr = isocorrmag - 1
    
    automag += zpt814
    sdss_r_auto = automag - 1
    
    print('Retrieving rmag from catalog')
    sdss_r_cat = ECO_phot_cat.rmag.loc[ECO_phot_cat.ECOID==obj].values[0]
    
    sdssr_petro_calc.append(sdss_r_petro)
    sdssr_iso_calc.append(sdss_r_iso)
    sdssr_isocorr_calc.append(sdss_r_isocorr)
    sdssr_auto_calc.append(sdss_r_auto)
    sdssr_cat.append(sdss_r_cat)
#    y_err.append(magerr)
    hdu_f814w_coadd.close()

os.chdir('..')
print('Plotting')
x = np.linspace(0,len(good_Obj_subset)+1,len(good_Obj_subset))
my_xticks = good_Obj_subset
fig1 = plt.figure(figsize=(10,8))
plt.xticks(x, my_xticks,rotation=90)
plt.scatter(x,sdssr_petro_calc, c='r',label='calculated petro rmag')
plt.scatter(x,sdssr_iso_calc, c='b',label='calculated iso rmag')
plt.scatter(x,sdssr_isocorr_calc, c='c',label='calculated iso corr rmag')
plt.scatter(x,sdssr_auto_calc, c='m',label='calculated auto rmag')
plt.scatter(x,sdssr_cat, c='g',label='catalog rmag')
plt.xlabel('ECOID')
plt.ylabel('rmag')
plt.gca().invert_yaxis()
plt.legend(loc='lower left')
plt.title('Comparison of calculated rmag and rmag from ECO photometric table')
plt.savefig('calcrmag_catrmag.png')
    
    
    