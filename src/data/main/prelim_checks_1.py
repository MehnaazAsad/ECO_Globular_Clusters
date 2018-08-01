# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 21:57:40 2017

@author: asadm2
"""
### DESCRIPTION
#This script takes in the original ECO catalog and does the following:
#adds a filter column, carries out an exposure time check, a "imaging in 
#at least 2 filters" check and produces 2 text files containing the list of 
#ECOIDs that passed these checks as well as the corresponding URLs for the
#images.

from cosmo_utils.utils import work_paths as cwpaths
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

### Paths
dict_of_paths = cwpaths.cookiecutter_paths()
path_to_raw = dict_of_paths['raw_dir']
path_to_interim = dict_of_paths['int_dir']
path_to_figures = dict_of_paths['plot_dir']
#path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
#'data/raw/'
#path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
#'data/interim/'

#Original ECO catalog
ECO = path_to_raw + 'Available_HST_Data_ECO.txt'
#Reading catalog
ECO = pd.read_csv(ECO, delimiter='\s+', header=None, names=['ECOID', 'HSTOBJ'\
      , 'RA', 'DEC', 'exptime', 'camera', 'filename'])

#Formating column values
ECO['exptime'] = pd.to_numeric(ECO['exptime'],errors='coerce')
ECO['filename'] = ECO['filename'].astype('str')

#Number of elements in ECO
files_arr = ECO['filename'].values
n_files = len(files_arr)

#To retrieve filter name from filename values and add filter column to existing 
#ECO catalog
wfc3_ir = ['f110w','f125w','f160w']   
wfc3_uvis = ['f606w','f600lp']
filters = [[] for x in range(n_files)]
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

filters_unique = np.unique(filters)
filters_counter = Counter(filters)

#Adding filter array to DataFrame
ECO.loc[:, 'filters'] = filters

#Creating filter and minimum exposure time dictionary for first check
#Calculated using HST ETC
exptime_arr = [9399, 3671, 3331, 1319, 2055, 2236, 1758, 10337, 2045, 1237, 
               2290, 3853, 1928101311, 73024829, 275363, 1241, 31705, 26575,
               6021, 3548, 3723, 2053, 2249, 3368, 5275, 4069, 171413, 31062,
               11431, 5789, 8520, 10071, 6677, 24445, 12605, 10757, 50294]
exp_fil_dict = dict(zip(filters_unique, exptime_arr ))


ECOID_groups = ECO.groupby(['filters','ECOID']) 
ECO_keys     = ECOID_groups.groups.keys()

#Here the for loop is going object by object and filter by filter and grabbing 
#all related images, calculating the sum of their individual exptime and 
#comparing that value to the minimum for that filter according to the
#dictionary
ECO_match = []
for key in ECO_keys:
    if ECOID_groups.get_group(key).exptime.sum() >= exp_fil_dict[key[0]]:
        ECO_match.append(key) #"good" keys
ECO_match = np.array(ECO_match)
ECO_unq_id = np.unique(ECO_match.T[1]) #transpose ECO_match and get the unique 
                                       #ECOIDs
#This tells us how many objects passed the exptime check
print('Unique ECO Objects: {0}'.format(ECO_unq_id.size))

### Total number of images that you can use

#get all matching rows and all corresponding columns
ECO_match_img = ECO.loc[ECO.ECOID.isin(ECO_unq_id), : ] 
#for each uniqID key get all images corresponding to that ID from main catalog 
print('Total images that you can use: {0}'.format(ECO_match_img.shape[0]))
###

### Total number of filters per ECOID

ECO_match_img_fil_arr = [[ecoid, ECO_match_img.filters.loc[ECO_match_img.\
                                                           ECOID==ecoid].size]\
    for ecoid in ECO_unq_id]
#For each unique ECOID get each filter image that matches it from the catalog 
#that passed the exptime cut
ECO_match_img_fil_arr = np.array(ECO_match_img_fil_arr)
#Next to each filter there is a number of how many images were taken using 
#that filter so if we sum the entire column we should get the total number of 
#images
ECO_match_img_fil_arr.T[1].astype(float).sum() 
###

### Number of images per ECOID
ECOmatch_ID_counter = Counter(ECO_match_img.ECOID) #also gives total number of 
                                                   #images
ID_freq = ECOmatch_ID_counter.most_common() #gives most common ECOIDs
fig1 = plt.figure(figsize=(10,8))
labels_ID, values = zip(*ID_freq)
indexes = np.arange(len(labels_ID))
plt.plot(indexes, values, 'ro-')
plt.title('Image count per object')
plt.xlabel('Object number') #get corresponding ECOIDs from labels_ID
plt.ylabel('Number of images')
plt.show()
###

### Number of images per filter
ECOmatch_filters_counter = Counter(ECO_match_img.filters)
Filters_freq = ECOmatch_filters_counter.most_common()

fig2 = plt.figure(figsize=(10,8))
labels_fil, values = zip(*Filters_freq)
indexes = np.arange(len(labels_fil))
my_xticks = np.array(labels_fil)
values = np.array(values)
plt.xticks(indexes, my_xticks, rotation=90)
plt.plot(indexes, values, 'ro-')
plt.title('Number of images taken using each filter')
plt.xlabel('Filter name')
plt.ylabel('Number of images')
plt.show()
###

#Checking how many objects have images taken using at least two filters
ECO_match2 = []
filtercount = [] 
counter = 0
obj_fil_arr = [[ecoid, ECO_match_img.filters.loc[ECO_match_img.ECOID==ecoid],]\
               for ecoid in ECO_unq_id]

obj_fil_arr_subset = []
for obj in obj_fil_arr:
    opt_fil_arr = []
    ir_fil_arr = []
    for fil in obj[1].values:
        for opt_fil in ['475','555','606']:
            if opt_fil in fil:
                opt_fil_arr.append(fil)        
        for ir_fil in ['815','850']:
            if ir_fil in fil:
                ir_fil_arr.append(fil)
    if len(opt_fil_arr) > 0 and len(ir_fil_arr) > 0:
        obj_fil_arr_subset.append((obj[0],opt_fil_arr+ir_fil_arr))

for i in range(len(obj_fil_arr)): #len(obj_fil_arr) same as ECO_unq_id.size
    filtercount.append(len(np.unique(obj_fil_arr[i][1])))
    if len(np.unique(obj_fil_arr[i][1])) >= 2:
        ECO_match2.append(obj_fil_arr[i][0])
        counter += obj_fil_arr[i][1].size #to check that len(ECO_match_img2) 
                                          #accounted for all images

### Histogram of number of filters each object has
fig3 = plt.figure(figsize=(10,8))
filtercount = np.array(filtercount)
plt.hist(filtercount, bins=np.arange(filtercount.min(), filtercount.max()+1))
plt.xlabel('Number of filters')
plt.title('Histogram of number of filters each object has')
### 
    
ECO_match2 = np.array(ECO_match2) #array of ECOIDs that passed second check
print('Total objects that have at least two filters: {0}'.format(ECO_match2.size))
print('Total number of images corresponding to 208 objects: {0}'.format(counter))

#Match and return all columns for this ECOID from ECO catalog
ECO_match_img2 = ECO.loc[ECO.ECOID.isin(ECO_match2),: ]
#final list of urls
ECO_match_url = np.array(ECO_match_img2.filename)  

#Array of ecoid with corresponding urls
obj_url_arr = [[ecoid, ECO_match_img2.filename.loc[ECO_match_img2.\
                                                   ECOID==ecoid],] \
    for ecoid in ECO_match2]

ID_arr = []
url_arr = []
for i in range(len(obj_url_arr)):
    ID_arr.append(obj_url_arr[i][0])
    url_arr.append(obj_url_arr[i][1])
    
ecoid_url_dict = dict(zip(ID_arr, url_arr)) #dict of IDs and urls        

### Text files

#Saving text file of ECOIDs
os.chdir(path_to_interim)
        
ecoidkeys = sorted(ecoid_url_dict.keys())
np.savetxt('Obj_arr.txt',ecoidkeys,fmt='%s')

#Saving text files of URLs for each ECOID in ecoidkeys
if not os.path.exists(path_to_interim + 'URLs'):
    os.mkdir('URLs')
os.chdir('URLs')
for i in range(len(ecoidkeys)):   
    np.savetxt(ecoidkeys[i] + '.txt',ecoid_url_dict[ecoidkeys[i]],fmt='%s')
###
    
'''
contents = pd.read_csv('ECO00649rev.txt',header=None)
contents.columns = ['filename']
#ECO2 = contents.filename.values

contents.filename = 'http://hla.stsci.edu/cgi-bin/' + contents.filename.astype(str)
ECO2 = ECO.loc[(ECO.filename.isin(contents.filename)) & (ECO.ECOID=='ECO00649'),: ]

ECOID_groups = ECO2.groupby('filters') #returns dict of grouped dataframes
ECO_keys     = ECOID_groups.groups.keys() #get the key list of each group in the dict
 
ECO_match3 = []
for key in ECO_keys:
    if ECOID_groups.get_group(key).exptime.sum() >= exp_fil_dict[key]: #exptime check
        ECO_match3.append(key) #"good" keys
ECO_match3 = np.array(ECO_match3)
if len(ECO_match3) >= 2:
    print ("{0} is good!".format(np.unique(ECO2.ECOID)))
else:
    print ("{0} is not good!".format(np.unique(ECO2.ECOID)))

#To check diff between initial Obj_arr and goodandbadObj because 2 objects were
#missing

goodObj = pd.read_csv('goodObj.txt',header=None)
badObj = pd.read_csv('badObj.txt',header=None)

filenames = ['goodObj.txt', 'badObj.txt']
with open('goodandbadObj.txt', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())

goodandbadObj = pd.read_csv('goodandbadObj.txt',header=None)
'''

