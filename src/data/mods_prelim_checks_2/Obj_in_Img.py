# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 21:13:49 2017

@author: asadm2
"""

### DESCRIPTION
#This script checks if the RA and DEC of the galaxy are in the image

from astropy.wcs import WCS
from astropy.io import fits
from pathlib2 import Path
from glob import glob
import os

path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'

def obj_in_img(RA,DEC,objname):
    
    """
    This function checks if the RA and DEC of the galaxy are in the image
    
    Args: 
        RA: RA of galaxy
        DEC: DEC of galaxy
        objname: ECOID of the galaxy
    
    Returns:
        Text file per object of all the filenames that passed this check
        Text file of errors from handling FITS file
    """
    arr = glob('*getdata*') #grab downloaded images
    
    for filename in arr:
        try:
            hdulist = fits.open(filename)
            header = hdulist[1].header
            w = WCS(header)
            px,py = w.wcs_world2pix(RA,DEC,1) #gets pixel position
            px = int(px)
            py = int(py)
            if px <= header['NAXIS1'] and px >= 0 and py <= header['NAXIS2']\
            and py >= 0:
                with open(objname+'_rev.txt', 'a') as newfile:
                    newfile.write(filename+'\n')
            else:
                os.remove(filename)
            hdulist.close()
        except KeyboardInterrupt:
            return
        #To catch empty or corrupt FITS files and other errors
        except Exception as e: 
            with open(path_to_interim+'Error_prelim2.txt', 'a') as newfile2:
                newfile2.write("{0} for object {1} raises an error: {2}".\
                               format(filename,objname,e))
                newfile2.write('\n')
    
    file_name = Path(objname+'_rev.txt')
    if file_name.is_file():
        newfile.close()
