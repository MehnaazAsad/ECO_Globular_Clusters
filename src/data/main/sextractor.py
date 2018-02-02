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

#good_Obj_subset = good_Obj[:4]['ECOID'].append(good_Obj[5:9]['ECOID']).values

good_Obj_subset = ['ECO00026'] #remove after testing

sdssr_stpetro_calc = []
sdssr_abpetro_calc = []
sdssr_iso_calc = []
sdssr_isocorr_calc = []
sdssr_auto_calc = []
sdssr_cat = []

#y_err = []
for index,obj in enumerate(good_Obj_subset):
    print('{0}/{1} {2}'.format(index+1,len(good_Obj_subset),obj))
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir('../../../data/raw/'+obj) #Change raw back to interim
    elif os.path.basename(dir_path) != obj:
        os.chdir('../'+obj)
        
#    comb_coadd = glob(obj+'_comb_coadd.fits')[0]
#    f814_coadd = glob(obj+'_acs_wfc_f814w_coadd.fits')[0]

    f814_coadd = 'hlsp_coma_hst_acs-wfc_v24_f814w_v1_ivm-drz-cl' #remove after testing

    hdu_f814w_coadd = fits.open(f814_coadd)
    
    data = hdu_f814w_coadd[1].data #remove 63-65 after testing
    hdr = hdu_f814w_coadd[1].header
    fits.writeto('test_'+f814_coadd, data, hdr,output_verify='ignore')
    
    
    print('Starting source extractor')
#    subprocess.call(['sex',comb_coadd+","+f814_coadd,'-ANALYSIS_THRESH','1.5',\
#    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
#    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1','-PIXEL_SCALE','0.0',\
#    '-CATALOG_NAME',obj+'_acs_wfc_f814w.cat']) #CHANGE THIS
    subprocess.call(['sex','test_'+f814_coadd,'-ANALYSIS_THRESH','1.5',\
    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1','-PIXEL_SCALE','0.0',\
    '-CATALOG_NAME',obj+'_acs_wfc_f814w.cat']) #CHANGE THIS
    print('Finished running source extractor') #remove 73-76 after testing
    
      
#    prihdr = hdu_f814w_coadd[0].header
    prihdr = hdu_f814w_coadd[1].header #delete after testing
    photflam814 = prihdr['PHOTFLAM']
    photzpt814 = prihdr['PHOTZPT']
    photplam814 = prihdr['PHOTPLAM']
    
    print('Calculating zeropoint')
    stzpt814 = (-2.5*np.log10(photflam814) + photzpt814)
    abzpt814 = -2.5*np.log10(photflam814) - (5*np.log10(photplam814)) - 2.408

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


#    cat_sextractor = SkyCoord(f814w_cat['x_world']*u.deg,\
#                              f814w_cat['y_world']*u.deg)
#    cat_eco = SkyCoord(ra*u.deg, dec*u.deg)
#    idx_sdss, d2d_sdss, d3d_sdss = cat_eco.match_to_catalog_sky(cat_sextractor)
#
#    f814mag = f814w_cat.petro_mag.values[idx_sdss]
#    isomag = f814w_cat.iso_mag.values[idx_sdss]
#    isocorrmag = f814w_cat.isocorr_mag.values[idx_sdss]
#    automag = f814w_cat.auto_mag.values[idx_sdss]
#    
#    petroflux = f814w_cat.petro_flux.values[idx_sdss]
    

    
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
    
    os.chdir('..')
    with open('sextractor_magflux_ECO00026.txt','a') as newfile: #change name
        newfile.write('{0},{1}\n'.format(f814mag,f814_coadd)) #remove name of image, add back petroflux
        newfile.close()
    os.chdir(obj)
    
    print(f814mag)

    print('Calculating rmag')
    f814stmag = f814mag + stzpt814
    f814abmag = f814mag + abzpt814
    sdss_r_petro_st = f814stmag - 0.9898
    sdss_r_petro_ab = f814abmag
    
#    isomag += stzpt814
#    sdss_r_iso = isomag - 1
#    
#    isocorrmag += stzpt814
#    sdss_r_isocorr = isocorrmag - 1
#    
#    automag += stzpt814
#    sdss_r_auto = automag - 1
    
    print('Retrieving rmag from catalog')
    sdss_r_cat = ECO_phot_cat.rmag.loc[ECO_phot_cat.ECOID==obj].values[0]
    
    sdssr_stpetro_calc.append(sdss_r_petro_st)
    sdssr_abpetro_calc.append(sdss_r_petro_ab)
    
#    sdssr_iso_calc.append(sdss_r_iso)
#    sdssr_isocorr_calc.append(sdss_r_isocorr)
#    sdssr_auto_calc.append(sdss_r_auto)
    sdssr_cat.append(sdss_r_cat)
#    y_err.append(magerr)
    hdu_f814w_coadd.close()

os.chdir('..')
print('Plotting')
x = np.linspace(0,len(good_Obj_subset)+1,len(good_Obj_subset))
my_xticks = good_Obj_subset
fig1 = plt.figure(figsize=(10,8))
plt.xticks(x, my_xticks,rotation=90)
plt.scatter(x,sdssr_stpetro_calc,s=50, c='r',label='calculated petro rmag (st)')
plt.scatter(x,sdssr_abpetro_calc,s=50, c='b',label='calculated petro rmag (ab)')
#plt.scatter(x,sdssr_iso_calc, c='b',label='calculated iso rmag')
#plt.scatter(x,sdssr_isocorr_calc, c='c',label='calculated iso corr rmag')
#plt.scatter(x,sdssr_auto_calc, c='m',label='calculated auto rmag')
plt.scatter(x,sdssr_cat,s=50, c='g',label='catalog rmag')
plt.xlabel('ECOID')
plt.ylabel('rmag')
plt.gca().invert_yaxis()
plt.legend(loc='lower left')
plt.title('Comparison of calculated rmag and rmag from ECO photometric table')
plt.savefig('calcrmag_catrmag_ECO00026.png') #change after testing
    
    
    