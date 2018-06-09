# -*- coding: utf-8 -*-
"""
@author: asadm2
"""
from glob import glob
import pandas as pd
import numpy as np
import shutil
import sys
import os

### DESCRIPTION
#This script downloads the images per object that passed prelim_checks_1 from 
#Obj_arr.txt. It then checks if the object is in the image andrepeats the 
#exposure time and "at least 2 filter" checks. The remaining "good" objects 
#that passed these checks are written to a text file and the "bad" objects 
#to a separate text file. 

### Importing mods
sys.path.insert(0, '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
                'data/mods_prelim_checks_2/')
from Obj_in_Img import obj_in_img
from Exp_fil2 import expfil2

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'

#If re-running script then these files have to be removed since
#they are being opened in 'append' mode below
if os.path.isfile(path_to_interim + 'Error_prelim2.txt'): 
    os.remove(path_to_interim + 'Error_prelim2.txt')
if os.path.isfile(path_to_interim + 'goodObj.txt'): 
    os.remove(path_to_interim + 'goodObj.txt')
if os.path.isfile(path_to_interim + 'badObj.txt'): 
    os.remove(path_to_interim + 'badObj.txt')

#Original ECO catalog
ECO = path_to_raw + 'Available_HST_Data_ECO.txt'
#Reading catalog
ECO = pd.read_csv(ECO, delimiter='\s+', header=None, \
                  names=['ECOID', 'HSTOBJ', 'RA', 'DEC', 'exptime', \
                         'camera', 'filename'])

#Formating column values
ECO['exptime'] = pd.to_numeric(ECO['exptime'],errors='coerce')
ECO['filename'] = ECO['filename'].astype('str')

#Reading text file of objects that passed prelim_checks_1
ECOIDs = pd.read_csv(path_to_interim + 'Obj_arr.txt', header=None,\
                     names=['ECO_ID'])
ECOID_arr = ECOIDs.ECO_ID.values

for obj in ECOID_arr:
    print("I am currently on {0}".format(obj))
    #Download images
    os.system("python pget_ECO.py ../../../data/interim/URLs/"+obj+'.txt')
    
    #Read RA and DEC from catalog
    RA = np.unique(ECO.RA.loc[ECO.ECOID==obj])[0]
    DEC = np.unique(ECO.DEC.loc[ECO.ECOID==obj])[0]
    obj_in_img(RA,DEC,obj) #Checks if galaxy is in the image
    dir_path = obj
    obj_list_path = obj+'_rev.txt' #Output from obj_in_img module
    
    if os.path.exists(obj_list_path):
        if not os.path.exists(dir_path):
            os.makedirs(dir_path) #make new folder where obj_rev.txt resides
        arr2 = glob('getdata*')
        
    for filename in arr2:
        #Move images to folder for object
        os.rename(filename, dir_path+"/"+filename)
        #Exposure time and "at least 2 filter" check
        Obj_success = expfil2(obj,obj_list_path) 
        
        if Obj_success:
            #Move obj_rev.txt file to folder for object
            os.rename(obj_list_path, dir_path+"/"+obj_list_path)      
            print('{0} is good'.format(obj))
            
        else:
            print('{0} is bad'.format(obj))
            #Delete both folder and obj_rev.txt file
            shutil.rmtree(obj)
            os.remove(obj_list_path)
    
    else:
        #If obj_rev.txt doesn't exist i.e. this ECOID didn't pass the "object
        #in image" check
        with open(path_to_interim + 'badObj.txt', 'a') as newfile:
            newfile.write(obj+'\n')
        print('{0} is bad'.format(obj)) 
        newfile.close()