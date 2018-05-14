# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:58:27 2017

@author: asadm2
"""
from __future__ import division, absolute_import, print_function
import pandas as pd
import numpy as np
import os
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from astropy.io import fits
import astropy
from glob import glob
from astropy.io.fits import getdata

def filter_append(filename):
    """
    This function retrieves the filter from the filename and appends it
    as a new column called 'filters' to the ECO catalog.
    
    Args: 
        filename: raw ECO catalog
    
    Returns:
        nothing but modifies raw ECO catalog   
    """
    ECO = pd.read_csv(filename, delimiter='\s+', header=None, names=['ECOID', 'HSTOBJ'\
          , 'RA', 'DEC', 'exptime', 'camera', 'filename'])
    
    ECO['exptime'] = pd.to_numeric(ECO['exptime'],errors='coerce')
    ECO['filename'] = ECO['filename'].astype('str')
    
    # Number of elements in ECO
    files_arr = ECO['filename'].values
    n_files = len(files_arr)
    
    wfc3_ir = ['f110w','f125w','f160w']   
    wfc3_uvis = ['f606w','f600lp']
    filters = [[] for x in xrange(n_files)]
    #Goes through every filename in catalog, retrieves filter name,
    #formats it and appends it to the filters array which is added 
    #in the end as a column to the ECO catalog.
    for i in range(len(ECO['filename'])):
        str_split = ECO['filename'][i].split(';')[1].split('_')
        filter_i = str_split[3]+'_'+str_split[4]+'_'+str_split[5]
        
        if 'ACS' in filter_i: #acs_wfc
            filter_i = filter_i.lower()
        
        elif 'd634' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
            
        elif 'm51' in filter_i: #acs-wfc fixed to acs_wfc       
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'tile' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'c_v' in filter_i: #acs-wfc fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            if 'acs-wfc' in filter_i:
                str_split = filter_i.split('-')
                filter_i = str_split[0]+'_'+str_split[1]
                
        elif 'ngc' in filter_i: #acs fixed to acs_wfc
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_wfc_'+str_split[2]
            
        elif '131009' in filter_i: #wfc3
            str_split = filter_i.split('_')
            if str_split[2] == 'f438w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            elif str_split[2] == 'f775w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
                                                 
        elif 'par' in filter_i:#and any(str in filter_i for str in wfc3_ir):
            str_split = filter_i.split('_')
            if str_split[2] == 'f606w':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            elif str_split[2] == 'f125w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f160w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f110w':
                filter_i = str_split[0]+'_ir_'+str_split[2]
            elif str_split[2] == 'f600lp':
                filter_i = str_split[0]+'_uvis_'+str_split[2]
            
        elif 'w_wf' in filter_i: #wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[1]
        
        elif 'lp_wf' in filter_i: #wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[1]
            
        elif 'n4496' in filter_i: #all wfpc2
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
            
            
        elif 'n5194' in filter_i: #all wfpc2 
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]
           
        elif 'u6614' in filter_i: #all wfpc2   
            str_split = filter_i.split('_')
            filter_i = str_split[0]+'_'+str_split[2]

        filters[i] = filter_i
        
    filters = np.asarray(filters)
    
    #Add the filters column to the raw ECO catalog
    ECO.loc[:, 'filters'] = filters

def make_dict(ECOnew,goodObj):
    """
    This function creates a dictionary of images grouped by their ECOID and
    filter used. 
    
    Args:
        ECOnew:  formatted ECO catalog with filter and new_filename column
        goodObj: list of ECO objects that meet the requirements for swarp-ing.
        
    Returns:
        The dictionary with two keys and the 'value' being the filename
        A list of keys sorted by ECOID
        The formatted ECO catalog
    """
    ECO_new  = pd.read_csv(ECOnew, delimiter='\t')
    
    good_Obj         = pd.read_csv(goodObj,header=None)
    good_Obj.columns = ['ECO_ID']
    #An array of the ECOIDs obtained from goodObj.txt
    arr1             = good_Obj.ECO_ID.values
    
    #Creates a new dataframe where you grab all the columns that match ECOID
    #in goodObj with ECOID in the formatted ECO catalog
    ECOnew_goodObj = ECO_new.loc[ECO_new.ECOID.isin(arr1), : ]
    
    #Group by ECOID and filters
    ECOnew_groups = ECOnew_goodObj.groupby(['ECOID','filters'])
    #Retrieve keys to the group sorted by ECOID
    ECOnew_keys   = sorted(ECOnew_groups.groups.keys())
    
    return ECOnew_groups, ECOnew_keys, ECO_new

def update_header(arr_imgs,obj1,filter_i):
    """
    This function takes an image and checks for CCDGAIN and EXPTIME keywords
    in the first 2 headers of the fits file. If it finds the keywords it 
    writes them to both headers if it is missing in either. It then writes
    the updated header to a new version of the fits file. After it has done 
    this for all images related to a filter, it writes the names of the files
    to a text file.
    
    Args:
        arr_imgs: list of images retrieved from dictionary using ECOID key
        obj1:     ECO object name
        filter_i: filter name for that ECO object 
        
    Returns:
        Either the text file with a list of names of the new version of fits 
        files
        OR
        an 'error' string if there weren't any images to write to the text
        file to begin with. 
    
    Raises:
        IOError: if 'empty or corrupt FITS file' or 'file does not exist'
    """
    for img in arr_imgs:
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        try:
            hdulist = fits.open(img,ignore_missing_end=True)
            #if there is only a primary header get the data from it
            if len(hdulist) == 1:
                data = getdata(img, 0, header=False)
            #if there is more than one header get data from the 'SCI' extension
            else:
                data = getdata(img, 1, header=False)
            #Get value of EXPTIME and PHOTZPT keyword from primary header and 
            #set CCDGAIN to a default value of 1
            EXPTIME = hdulist[0].header['EXPTIME']
            PHOTZPT = hdulist[1].header['PHOTZPT']
            CCDGAIN = 1.0
            #First pass locating value for gain
            for i in range(2):
                if len(hdulist) == 1:
                    break
                #Go through primary and secondary header and ignore the 
                #BinTable formatted header
                if not isinstance(hdulist[i],astropy.io.fits.hdu.table.BinTableHDU):
                    if 'CCDGAIN' in hdulist[i].header:
                        CCDGAIN = hdulist[i].header['CCDGAIN']
                        break
                    if 'GAIN' in hdulist[i].header:
                        CCDGAIN = hdulist[i].header['GAIN']
                        break
                    if 'ATODGAIN' in hdulist[i].header:
                        CCDGAIN = hdulist[i].header['ATODGAIN']
                        break
            
            print('Doing BUNIT check')
            for i in range(2):
                if len(hdulist) == 1:
                    bunit = hdulist[0].header['D001OUUN']
                    print('BUNIT was {0}'.format(bunit))
                    if bunit == 'counts':
                        data = data/EXPTIME
                        ZPT_NEW = 30
                        pixmod = 10**(-0.4*(PHOTZPT-ZPT_NEW))
                        data = data*pixmod
                        hdulist[0].header.set('BUNIT','COUNTS/S')
                        hdulist[0].header.set('MAGZPT',ZPT_NEW)
                        print('BUNIT is {0}'.format(hdulist[0].header['BUNIT']))
                else:
                    if 'BUNIT' in hdulist[i].header:
                        bunit = hdulist[i].header['BUNIT']
                        print('BUNIT was {0}'.format(bunit))
                        print('PHOTZPT was {0}'.format(PHOTZPT))
                        if bunit == 'COUNTS':
                            data = data/EXPTIME
                        if bunit == 'ELECTRONS':
                            data = data/(CCDGAIN*EXPTIME)
                            ZPT_NEW = 30
                            pixmod = 10**(-0.4*(PHOTZPT-ZPT_NEW))
                            data = data*pixmod
                        if bunit == 'ELECTRONS/S':
                            data = data/CCDGAIN
                            ZPT_NEW = 30
                            pixmod = 10**(-0.4*(PHOTZPT-ZPT_NEW))
                            data = data*pixmod
                        if bunit == 'ELECTRONS/SEC':
                            data = data/CCDGAIN
                        hdulist[i].header['BUNIT'] = 'COUNTS/S'
                        hdulist[i].header.set('MAGZPT',ZPT_NEW)
                        print('BUNIT is {0}'.format(hdulist[i].header['BUNIT']))
                        print('PHOTZPT is {0}'.format(hdulist[i].header['MAGZPT']))
            print('Done changing BUNIT')
            
            #Second pass to assign gain and exptime to headers
            for i in range(2):
                if len(hdulist) == 1:
                    break
                if not isinstance(hdulist[i],astropy.io.fits.hdu.table.BinTableHDU):
                    if 'CCDGAIN' not in hdulist[i].header:
                        hdulist[i].header.set('CCDGAIN',CCDGAIN)
                    if 'EXPTIME' not in hdulist[i].header:
                        hdulist[i].header.set('EXPTIME',EXPTIME)
            #If there was only one header write that header's data to new
            #version of fits image
            
            #Make new versions in interim/obj1 folder
            os.chdir('../../interim/'+obj1)
            print(img)
            img = os.path.splitext(img)[0]
            print(img)
            if len(hdulist) == 1:
                fits.writeto(img+'_test_'+filter_i+'.fits',data,hdulist[0].header,output_verify='ignore')
            #Else write the 'SCI' header's data to new version of fits image
            else:
                fits.writeto(img+'_test_'+filter_i+'.fits',data,hdulist[1].header,output_verify='ignore')
            hdulist.close()
            os.chdir('../../raw/'+obj1)
        #This is to catch 'empty or corrupt FITS file' or any other IOError
        #and write it to a text file along with the object name and the 
        #filter name
        except IOError as e:
            os.chdir('..')
            dir_path = os.getcwd()
            if os.path.basename(dir_path) == 'raw':
                os.chdir('../interim')
            with open('Error_swarpfil.txt','a') as newfile:              
                newfile.write('Object {0} and image {1} raises {2} error'.format\
                (obj1,img,e))
                newfile.write('\n')
                newfile.close()
            os.chdir('../raw/'+obj1)
    #For this object and filter combination grab all the new versions made
    os.chdir('../../interim/'+obj1)
    arr = glob('*test_'+filter_i+'.fits')
    print(len(arr))
    if len(arr) >= 1: #avoid empty cases where files have been removed earlier
                      #or don't exist at all since the dictionary also contains
                      #pairs of objects and filters that didn't meet the swarp
                      #requirements (didn't pass prelimary exptime or filter
                      #checks so those folders/images don't exist
                      
        #If new versions exist then write their names to a text file 
        with open(filter_i+'_img_list_testfil.txt','wb') as newfile2:
            for obj in arr:
                newfile2.write(obj)
                newfile2.write('\n')
            newfile2.close()
        #If text file exists return the name
        return filter_i+'_img_list_testfil.txt'
    #If text file doesn't exist return this string
    return 'error'

def percent_blank(coadd,obj,filter_i):
    """
    This function takes in a coadd, gets an array of all pixel values of the
    image and calculates the percentage of the image that is blank i.e. any 
    pixel that has a value of exactly 0
    
    Args: 
        coadd:    coadded fits image made using swarp from images from the same 
                  filter
        obj:      ECO object name
        filter_i: filter name for that ECO object
    
    Returns:
        nothing but creates a new text file
        
    """
    #Get array of pixel values
    scidata = coadd[0].data 
    #Flatten n-d array to 1-d 
    new  = np.ravel(scidata)
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
    #Move back up to interim folder
    os.chdir('..')
    #Write the percentage blank along with the object name and filter 
    #name to a text file
    with open('percent_blankv2.txt','a') as newfile:
        newfile.write('{0},{1},{2}%\n'.format(obj,filter_i,percentage))
        newfile.close()

    

        
    