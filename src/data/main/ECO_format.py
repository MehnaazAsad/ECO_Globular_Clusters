# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 22:21:18 2017

@author: asadm2
"""
### DESCRIPTION
#This script reformats the original ECO catalog to add a column for the 
#filter that the image was taken using as well as a column for a new version
#of the filename in the catalog (without http...=). This script returns a 
#the newly formatted catalog as a new text file (ECO_formatted.txt) which will
#be used from now on instead of the original ECO catalog. Only formatting 
#changes have been made.

import numpy as np
import imp

### Importing module
foo = imp.load_source('filter_append', '../mods_ECO_format/Filter_append.py')

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'

### Original ECO catalog
ECO = path_to_raw + 'Available_HST_Data_ECO.txt'

ECO = foo.filter_append(ECO)
arr = ECO.filename.values

### Changing filenames of images in the catalog
new_filename = []
for filename in arr:
    new_name = filename.split('=')[2]
    new_filename.append(new_name)

new_filename = np.asarray(new_filename)
ECO.loc[:, 'new_filename'] = new_filename

ECO.to_csv(path_to_interim+'ECO_formatted.txt',index=False, sep='\t',\
           columns=['ECOID', 'HSTOBJ', 'RA', 'DEC', 'exptime', 'camera', \
                    'filename','filters','new_filename'])
