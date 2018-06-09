#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 23:25:33 2017

@author: asadm2
"""
### DESCRIPTION
#This script runs SExtractor, matches the object with its detection in the SE
#catalog and retrieves the petrosian magnitude. Magnitude conversions are 
#applied (NEED TO BE CHANGED SINCE ZEROPOINT WAS RESCALED) and comparisons
#are made with the sdss magnitude from the ECO photometry catalog. This script
#returns 3 text files that are described below.

from __future__ import division, print_function, absolute_import
from astropy.coordinates import SkyCoord
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import rc
from glob import glob
import pandas as pd
import numpy as np
import subprocess
import os

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'
path_to_figures = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'reports/figures/'

##Formatting setup for plots
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True,fontsize=16)
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#If re-running script then these files hav to be removed since
#they are being opened in 'append' mode in this script
if os.path.isfile(path_to_interim + 'magnitude_errors_degmatch.txt'): 
    os.remove(path_to_interim + 'magnitude_errors_degmatch.txt')
if os.path.isfile(path_to_interim + 'catmatch_separation.txt'): 
    os.remove(path_to_interim + 'catmatch_separation.txt')
if os.path.isfile(path_to_interim + 'magnitude_errors_pixelmatch.txt'): 
    os.remove(path_to_interim + 'magnitude_errors_pixelmatch.txt')


#Comacandidates.txt is a file of all ECOIDs that fall in the halo group of mass
#10**14.5 solar masses from Figure 1 in ecodr1.py
comaObj = path_to_interim + 'comacandidates.txt'
coma_Obj = pd.read_csv(comaObj,header=None, names = ['ECOID'])
arr_comaObj = coma_Obj.ECOID.values

#ECO catalog
ECOnew  = path_to_interim + 'ECO_formatted.txt'
ECO_new  = pd.read_csv(ECOnew, delimiter='\t')

#ECO photometry catalog
ECOphotcat = path_to_interim + 'ECOphot_final.txt'
ECO_phot_cat  = pd.read_csv(ECOphotcat, delim_whitespace=True,header=None,\
                            skiprows=1,names=['ECOID','RA','DEC','VLG','umag',\
                                              'gmag','rmag','imag','zmag',\
                                              'Jmag','Hmag','Kmag'])

#subset of coma candidates I was testing that had imaging in f814w and f475w
coma_Obj_subset = coma_Obj[:4]['ECOID'].append(coma_Obj[5:9]['ECOID']).values

#Initializing arrays
sdssr_stpetro_calc = []
sdssr_abpetro_calc = []
sdssi_stpetro_calc = []
sdssr_cat = []
sdssi_cat = []

objs_to_plot = []
d2d_arr = []

for index,obj in enumerate(coma_Obj_subset): #change back to arr_comaObj once
                                             #done testing this subset
    print('{0}/{1} Object: {2}'.format(index+1,len(arr_comaObj),obj))
    dir_path = os.getcwd()
    if os.path.basename(dir_path) == 'main':
        os.chdir(path_to_interim + obj)
    elif os.path.basename(dir_path) != obj:
        os.chdir('../'+obj)
    
    #Grab combined coadd, and coadds in f814w and f475w filters per object
    comb_coadd = glob(obj+'_comb_coadd.fits')[0]
    f814_coadd = glob(obj+'_acs_wfc_f814w_coadd.fits')[0]
    f475_coadd = glob(obj+'_acs_wfc_f475w_coadd.fits')[0]
    
    #Open images
    hdu_f814w_coadd = fits.open(f814_coadd)
    hdu_f475w_coadd = fits.open(f475_coadd)
    
#    data = hdu_f814w_coadd[1].data #remove 63-65 after testing
#    hdr = hdu_f814w_coadd[1].header
#    fits.writeto('test_'+f814_coadd, data, hdr,output_verify='ignore')
    
    #Begin SExtractor
    #It outputs a catalog with whatever name you pass as the value of CATALOG_
    #NAME parameter
    print('Starting source extractor')
    subprocess.call(['sex',comb_coadd+","+f814_coadd,'-ANALYSIS_THRESH','1.5',\
    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1','-PIXEL_SCALE','0.0',\
    '-CATALOG_NAME',obj+'_acs_wfc_f814w.cat'])
    subprocess.call(['sex',comb_coadd+","+f475_coadd,'-ANALYSIS_THRESH','1.5',\
    '-BACK_SIZE','128','-DEBLEND_MINCONT','0.0025', '-DETECT_THRESH','1.5',\
    '-DETECT_MINAREA','9','-SEEING_FWHM','0.1','-PIXEL_SCALE','0.0',\
    '-CATALOG_NAME',obj+'_acs_wfc_f475w.cat'])
    
    #Retrieve PHOTFLAM, PHOTZPT and PHOTPLAM keywords from headers of both
    #f814w and f475w coadds
    prihdrf814 = hdu_f814w_coadd[0].header
    photflam814 = prihdrf814['PHOTFLAM']
    photzpt814 = prihdrf814['PHOTZPT']
    photplam814 = prihdrf814['PHOTPLAM']
    
    prihdrf475 = hdu_f475w_coadd[0].header
    photflam475 = prihdrf475['PHOTFLAM']
    photzpt475 = prihdrf475['PHOTZPT']
    photplam475 = prihdrf475['PHOTPLAM']
    
    #OLD WAY OF CALCULATION BEFORE ZEROPOINT WAS SCALED. NOW YOU JUST NEED TO 
    #ADD 30 TO WHATEVER PETRO MAG SEXTRACTOR RETURNS
    print('Calculating zeropoint')
    ###REPLACE FROM HERE
    stzpt814 = (-2.5*np.log10(photflam814) + photzpt814)
    abzpt814 = -2.5*np.log10(photflam814) - (5*np.log10(photplam814)) - 2.408
    
    stzpt475 = (-2.5*np.log10(photflam475) + photzpt475)
    abzpt475 = -2.5*np.log10(photflam475) - (5*np.log10(photplam475)) - 2.408
    ###TO HERE
    
    #Get RA and DEC of object
    ra = np.unique(ECO_new.RA.loc[ECO_new.ECOID==obj])[0]
    dec = np.unique(ECO_new.DEC.loc[ECO_new.ECOID==obj])[0]
    
    #Converting to pixel coordinates
    print('Converting to pixel coordinates')
    wcs = WCS(prihdrf814)   
    print(wcs)
    xx,yy = wcs.wcs_world2pix(ra,dec,1)
    print(xx,yy)

    print('Reading SE catalog')
    
    #Reading SE catalogs from both filters
    
    #'If-else' required because for some objects the parameters asked for, as 
    #specified in default.param file, were more than others. PROBABLY BEST 
    #TO MAKE THIS MORE UNIFORM AND HAVE THE SAME PARAMETERS SELECTED FOR ALL
    #OBJECTS IN DEFAULT.PARAM FILE
    if obj in arr_comaObj.ECOID[:8].values:
        f814w_cat = pd.read_csv(obj+'_acs_wfc_f814w.cat',header=None,\
                                delim_whitespace=True,\
                                names=['iso_mag','isocorr_mag','auto_mag',\
                                       'petro_flux','petro_fluxerr',\
                                       'petro_mag','petro_magerr',\
                                       'petro_radius','xmin_image',\
                                       'ymin_image','xmax_image','ymax_image',\
                                       'x_image','y_image','x_world',\
                                       'y_world','a_image','class_star'],\
                                       comment='#')
                                       
        f475w_cat = pd.read_csv(obj+'_acs_wfc_f475w.cat',header=None,\
                                delim_whitespace=True,\
                                names=['iso_mag','isocorr_mag','auto_mag',\
                                       'petro_flux','petro_fluxerr',\
                                       'petro_mag','petro_magerr',\
                                       'petro_radius','xmin_image',\
                                       'ymin_image','xmax_image','ymax_image',\
                                       'x_image','y_image','x_world',\
                                       'y_world','a_image','class_star'],\
                                       comment='#')
    else:
        f814w_cat = pd.read_csv(obj+'_acs_wfc_f814w.cat',header=None,\
                                delim_whitespace=True,\
                                names=['petro_mag',\
                                       'petro_magerr','petro_radius',\
                                       'xmin_image','ymin_image','xmax_image',\
                                       'ymax_image','x_image','y_image',\
                                       'x_world','y_world','a_image',\
                                       'class_star'],comment='#')        
    
    #Retrieving petro mag for object from catalog
    #Two ways: 1. Use SkyCoord's match_to_catalog_sky to find closest RA and 
    #DEC match between actual object and those detected by SE. 2. Use pixel 
    #coordinates of object and find which SE aperture the coordinates are in
    #by using pixel max and min coordinates in x and y for detected objects.
    print('Getting petro mag from SE catalog')
    
    ###METHOD 1
    try:
        #Select position columns in both catalogs and specify units (deg)
        cat_sextractor = SkyCoord(f814w_cat['x_world']*u.deg,\
                                  f814w_cat['y_world']*u.deg)
        cat_eco = SkyCoord(ra*u.deg, dec*u.deg)
        #Match both catalogs. idx_se are the indecies into cat_sextractor that 
        #get the closest matches, while d2d and d3d are the on-sky and 
        #real-space distances between the matches.
        idx_se, d2d_se, d3d_se = cat_eco.match_to_catalog_sky(cat_sextractor)
        d2d_arr.append(d2d_se)
        #Get petrosian magnitude from catalog using idx
        f814mag = f814w_cat.petro_mag.values[idx_se]
        
        #For some detections petro mag in SE catalog was 99.0
        if f814mag == 99.0:
            raise ValueError
            
    #Record this ValueError in text file    
    except ValueError as valueerror:
        os.chdir(path_to_interim)
        with open('magnitude_errors_degmatch.txt','a') as newfile:
            newfile.write('{0} has a mag of 99\n'.format(obj))
            newfile.close()
        print('Magnitude is 99. Moving to next object.')
        os.chdir(obj)
        hdu_f814w_coadd.close()
        continue
    
    #Recording separation between matches in text file
    print('Writing separation to file')
    os.chdir(path_to_interim)
    with open('catmatch_separation.txt','a') as newfile:
        newfile.write('Separation for {0} is: {1}'.format(obj,\
                      d2d_se.arcsec[0]))
        newfile.write('\n')
        newfile.close()
    os.chdir(obj)
    ###END OF METHOD 1   
    
    ###METHOD 2
    try:
        #Check if pixel position of object falls within any of the SE
        #detected apertures using their pixel ranges
        f814mag = f814w_cat.petro_mag.loc[((f814w_cat.xmin_image < [xx])&([xx]\
                                           < f814w_cat.xmax_image))&\
                                          ((f814w_cat.ymin_image < [yy])&([yy]\
                                           < f814w_cat.ymax_image))].values[0]
                                          
        f475mag = f475w_cat.petro_mag.loc[((f475w_cat.xmin_image < [xx])&([xx]\
                                           < f475w_cat.xmax_image))&\
                                          ((f475w_cat.ymin_image < [yy])&([yy]\
                                           < f475w_cat.ymax_image))].values[0]                                               
        
        #For some detections petro mag in SE catalog was 99.0
        if f814mag == 99.0:
            raise ValueError
        if f475mag == 99.0:
            raise ValueError
            
    #Record IndexError in text file
    except IndexError as indexerror:
        os.chdir(path_to_interim)
        with open('magnitude_errors_pixelmatch.txt','a') as newfile:
            newfile.write('{0} gives {1}\n'.format(obj,indexerror))
            newfile.close()
        print('Error in getting magnitude. Moving to next object.')
        os.chdir(obj)
        hdu_f814w_coadd.close()
        hdu_f475w_coadd.close()
        continue
    
    #Record ValueError in text file    
    except ValueError as valueerror:
        os.chdir(path_to_interim)
        with open('magnitude_errors_pixelmatch.txt','a') as newfile:
            newfile.write('{0} has a mag of 99\n'.format(obj))
            newfile.close()
        print('Magnitude is 99. Moving to next object.')
        os.chdir(obj) 
        hdu_f814w_coadd.close()
        hdu_f475w_coadd.close()
        continue
    ###END OF METHOD 2
    
    #Retrieving errors in magnitudes                                         
    f814w_magerr = f814w_cat.petro_magerr.loc[((f814w_cat.xmin_image < \
                                                [xx])&([xx] < \
                                                f814w_cat.xmax_image))&\
                                              ((f814w_cat.ymin_image < \
                                                [yy])&([yy] < \
                                                f814w_cat.ymax_image))]\
                                                .values[0]

    f475w_magerr = f475w_cat.petro_magerr.loc[((f475w_cat.xmin_image < \
                                                [xx])&([xx] < \
                                                f475w_cat.xmax_image))&\
                                              ((f475w_cat.ymin_image < \
                                                [yy])&([yy] < \
                                                f475w_cat.ymax_image))]\
                                                .values[0]   
    
    f814mag = pd.to_numeric(f814mag)
    f475mag = pd.to_numeric(f475mag)
        
    print('f814w magnitude: {0}'.format(f814mag))
    print('f475w magnitude: {0}'.format(f475mag))
    
    #NOT VALID ANYMORE SINCE RESCALING OF ZEROPOINT
    print('Calculating rmag')
    ###MAY HAVE TO CHANGE FROM HERE
    f814stmag = f814mag + stzpt814
    f475stmag = f475mag + stzpt475
    
    f814abmag = f814mag + abzpt814
    
    #Using Lupton et al. transformations. Requires double checking!
    #Calculating apparent magnitudes to compare with those from the photometry 
    #catalog of ECO depending on which sdss filters f475w and f814w are closest
    #to.
    sdss_r_petro_st = f814stmag - 0.489 #Using average r-i for ECO of 0.289
    sdss_r_petro_ab = f814abmag
    sdss_i_petro_st = f814stmag - 0.748 #Using average i-z for ECO of 0.226
    ###TO HERE
    
    #Retrieving magnitudes of object from ECO photometry catalog
    print('Retrieving rmag from catalog')
    sdss_r_cat = ECO_phot_cat.rmag.loc[ECO_phot_cat.ECOID==obj].values[0]
    sdss_i_cat = ECO_phot_cat.imag.loc[ECO_phot_cat.ECOID==obj].values[0]
    
    #Appending to arrays initialized above
    sdssr_stpetro_calc.append(sdss_r_petro_st)
    sdssr_abpetro_calc.append(sdss_r_petro_ab)
    sdssi_stpetro_calc.append(sdss_i_petro_st)   
    sdssr_cat.append(sdss_r_cat)
    sdssi_cat.append(sdss_i_cat)
    objs_to_plot.append(obj)

    hdu_f814w_coadd.close()
    hdu_f475w_coadd.close()


###FIGURES
os.chdir(path_to_figures)
print('Plotting')

#Figure 1: comparison of SE retrieved r band magnitude and that from photometry 
#catalog
x = np.linspace(0,len(objs_to_plot)+1,len(objs_to_plot))
my_xticks = objs_to_plot

fig1 = plt.figure(figsize=(10,8))
plt.xticks(x, my_xticks,rotation=90)
for label,x,y in zip(objs_to_plot,sdssr_cat,sdssr_stpetro_calc):
    if label=='ECO00026':
        plt.scatter(x,y,s=50,c='red',label='ECO00026')
    else:
        plt.scatter(x,y,s=50,c='lightgrey')

#Plotting
plt.scatter(sdssr_cat,sdssr_stpetro_calc,s=50,c='lightgrey',\
            label='r mag comparison')
#Plotting a 1-1 line
plt.plot(sdssr_cat,sdssr_cat,'k--')
# =============================================================================
# ##Uncomment lines below if you want to compare with other magnitudes that
# ##SE returns or magnitudes in the ab system.
# plt.scatter(x,sdssr_abpetro_calc,s=50, c='b',\
# label='calculated petro rmag (ab)')
# plt.scatter(x,sdssr_iso_calc, c='b',label='calculated iso rmag')
# plt.scatter(x,sdssr_isocorr_calc, c='c',label='calculated iso corr rmag')
# plt.scatter(x,sdssr_auto_calc, c='m',label='calculated auto rmag')
# plt.scatter(x,sdssr_cat,s=50, c='g',label='catalog rmag')
# =============================================================================
plt.xlabel(r'ECO sdss mag')
plt.ylabel(r'SE sdss mag')
plt.legend(loc='best')
plt.title(r'Comparison of calculated and catalogued sdss r magnitudes',\
          fontsize=16)
plt.savefig('calcrmag_catrmag_coma.png')

#Figure 2: comparison of SE retrieved i band magnitude and that from photometry 
#catalog 
fig2 = plt.figure(figsize=(10,8))    
plt.scatter(sdssi_cat,sdssi_stpetro_calc,s=50,c='lightgrey',\
            label='i mag comparison')
plt.plot(sdssi_cat,sdssi_cat,'k--')
plt.xlabel(r'ECO sdss mag')
plt.ylabel(r'SE sdss mag')
plt.legend(loc='best')
plt.title(r'Comparison of calculated and catalogued sdss i magnitudes',\
          fontsize=16)
plt.savefig('calcimag_catimag_coma.png')

#Figure 3: Histogram of separations for METHOD 1
fig3 = plt.figure(figsize=(10,8))
plt.hist(d2d_arr, histtype='step')
plt.xlabel(r'separation (degrees)')
plt.savefig('catmatch_separation.png')  