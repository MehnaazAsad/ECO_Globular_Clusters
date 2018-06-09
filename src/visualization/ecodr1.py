#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 04:27:25 2017

@author: asadm2
"""

###DESCRIPTION
#This script plots current sample of objects out of the entire ECO catalog,
#separates this plot into single_halos (group mass of less than 10**14 solar
#masses) and coma_halos (group mass of more than 10**14 solar masses). Also,  
#plots velocity dispersion vs projected distance of the coma_halos and 
#expected GC count based on parent galaxy stellar mass for single_halos. 
#Note: The expected GC count is also calculated for coma_halos but not plotted. 

from astropy.coordinates import SkyCoord
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
from astropy import units as u
from matplotlib import rc
import pandas as pd
import numpy as np

### Paths
path_to_raw = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/raw/'
path_to_interim = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'data/interim/'
path_to_figures = '/fs1/masad/Research/Repositories/ECO_Globular_Clusters/'\
'reports/figures/'

#Formatting setup for plots
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True,fontsize=25)
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

#Reading original ECO catalog, photometry catalog, and most updated array of 
#'good' objects 
eco_dr1 = path_to_interim + 'eco_dr1.txt'
eco_phot = path_to_interim + 'ECOphot_final.txt'
goodObj = path_to_interim + 'goodObjv2.txt'

ECO_original  = pd.read_csv(eco_dr1, delimiter='\s+',header=None,usecols=\
                            [0,1,2,3,4,5,13,15],\
                            names=['ECOID','RA','DEC','velocity','M_r',\
                                   'log_stellarmass','groupcz',\
                                   'log_halogroupmass']) 

ECO_photometric = pd.read_csv(eco_phot, delimiter='\s+',skiprows=1,\
                              header=None,names=\
                              ['ECOID','RA','DEC','VLG','umag','gmag','rmag',\
                               'imag','zmag','Jmag','Hmag','Kmag'])

good_Obj = pd.read_csv(goodObj,header=None, names = ['ECOID'])


#Matching ECOIDs in good_Obj with those from the original ECO catalog and 
#extracting all columns from ECO catalog 

#NOTE! 
#2 ECOIDs are missing that were present in good_Obj: ECO12932 and ECO13577
ECO_goodObj_match  = ECO_original.loc[(ECO_original.ECOID.isin\
                                       (good_Obj.ECOID)),: ]

#Using equation 1 from Zaritsky et al. 2015 to calculate expected number of 
#GCs given stellar mass of galaxy
T_N = []
for val in ECO_goodObj_match.log_stellarmass:
    val = 10**val
    if 10**8.5 <= val < 10**10.5:
        T_n = (10**6.7)*(val**(-0.56))
    elif 10**10.5 <= val <= 10**11:
        T_n = 8.3
    elif val > 10**11:
        T_n = (10**(-6.11))*(val**(0.63))
    T_N.append(T_n)

#Divide out by 10**9 to go from specific frequency to frequency    
T_N = np.log10(T_N*((10**ECO_goodObj_match.log_stellarmass)/10**9))

#Figure 1: Current sample (red) over entire ECO catalog (grey)
fig1 = plt.figure(figsize=(10,8))
plt.scatter(ECO_original.log_halogroupmass,ECO_original.M_r, c='lightgrey',\
            label='catalog')
plt.scatter(ECO_goodObj_match.log_halogroupmass,ECO_goodObj_match.M_r, c='r',\
            label='sample')
plt.xlabel(r'$\log\ M_h/M_\odot$')
plt.ylabel(r'$M_r$')
plt.gca().invert_yaxis()
plt.legend()
plt.title(r'Sample of ECO',fontsize=25)

#These are the red points from Figure 1 that are in groups of mass less than 
#10**14 solar masses
single_halos = ECO_goodObj_match.loc[ECO_goodObj_match.log_halogroupmass <14]

#Using equation 1 from Zaritsky et al. 2015 to calculate expected number of 
#GCs given stellar mass of galaxy
T_N_singlehalos = []
T_N_singlehalos_name = []
for val in single_halos.log_stellarmass:
    name = single_halos.ECOID.loc[single_halos.log_stellarmass==val].values[0]
    val = 10**val
    if 10**8.5 <= val < 10**10.5:
        T_n = (10**6.7)*(val**(-0.56))
    elif 10**10.5 <= val <= 10**11:
        T_n = 8.3
    elif val > 10**11:
        T_n = (10**(-6.11))*(val**(0.63))
        
    T_N_singlehalos.append(T_n)
    T_N_singlehalos_name.append(name)
 
#Divide out by 10**9 to go from specific frequency to frequency    
T_N_singlehalos = np.log10(T_N_singlehalos*((10**single_halos.log_stellarmass)\
                                            /10**9))

GC_count = pd.DataFrame({'ECOID': T_N_singlehalos_name, 'log_GC_count': \
                     T_N_singlehalos.values})

#Figure 2: Plot of number of GCs vs stellar mass of parent galaxy
fig2 = plt.figure(figsize=(10,8))
plt.scatter(single_halos.log_stellarmass,T_N_singlehalos,\
            c = single_halos.log_halogroupmass,cmap='plasma',s=100)
cbar = plt.colorbar()
cbar.set_label('Halo mass', rotation=270,labelpad=20,fontsize=20)
plt.xlabel(r'$\log\ M_*/M_\odot$')
plt.ylabel(r'$\log\ N$')
plt.title(r'Number of GCs vs parent galaxy stellar mass',fontsize=20)
plt.show()
#==============================================================================
#Figure 3: velocity dispersion vs projected distance

#These are the red points from Figure 1 that are in groups of mass greater than 
#10**14 solar masses thought to be part of the Coma cluster
coma_halos = ECO_goodObj_match.loc[ECO_goodObj_match.log_halogroupmass > 14]

CZ_coma = 6925 #km/s
Distance_coma = 102.975/0.705 #distance to Coma in Mpc/h
deltav = coma_halos.velocity - CZ_coma

#Calculating angular separation between Coma and coma_halos in order 
#to calculate projected distance
RADEC_coma = SkyCoord('12h59m48.7s', '+27d58m50s')
RADEC_eco_deg = SkyCoord(ra=coma_halos.RA.values*u.degree, dec=coma_halos.\
                        DEC.values*u.degree)
RADEC_eco = RADEC_eco_deg.to_string('hmsdms')
angular_separations_arr = []
for value in RADEC_eco:
    RA = value.split( )[0]
    DEC = value.split( )[1]
    value = SkyCoord(RA,DEC)
    angular_separation = RADEC_coma.separation(value)
    angular_separation = angular_separation.radian
    angular_separations_arr.append(angular_separation)
    
angular_separations_arr = np.array(angular_separations_arr)

projected_distance = (angular_separations_arr*Distance_coma) #in Mpc

#Axis labels and point labels (if necessary)
xlabel = r'\boldmath$r_{p} \left[h^{-1} Mpc \right]$'
ylabel = r'\boldmath$\Delta v\ \left[km/s\right]$'
labels = []
for name in coma_halos.ECOID.values:
    name = name.strip('ECO')
    labels.append(name)
labels = np.array(labels)

#g-r colour of objects
g_r = []
for name in coma_halos.ECOID.values:
    gmag = ECO_photometric.gmag.loc[ECO_photometric.ECOID==name].values[0]
    rmag = ECO_photometric.rmag.loc[ECO_photometric.ECOID==name].values[0]
    colour = gmag-rmag
    g_r.append(colour)
    

fig3 = plt.figure(figsize=(10,8))
plt.scatter(projected_distance,deltav,s=70,c=g_r)
cbar = plt.colorbar()
cbar.set_label('g-r', rotation=270,labelpad=20,fontsize=12)
plt.clim(0,1)
#for i, txt in enumerate(labels):
#    plt.annotate(txt, (projected_distance[i],deltav.values[i]),fontsize=5)
plt.axhline(y=1000,linestyle='dashed',color='r') #velocity dispersion of coma
plt.axhline(y=-1000,linestyle='dashed',color='r')
plt.axvline(x=2,linestyle='dashed',color='r')
plt.xlim(-0.01,)
plt.text(1.5, 650, r'$\mathrm{\frac{N_{coma}}{total}=\frac{30}{42}}$',\
         fontsize=25)
## Shading
x_arr = np.arange(-0.01,2.01, 0.1)
y1_arr = np.ones(len(x_arr ))*1000
y2_arr = np.ones(len(x_arr ))*-1000
plt.fill_between(x_arr, y1_arr, y2_arr, color='grey', alpha=0.2)
##

plt.xlabel(xlabel, fontsize=20)
plt.ylabel(ylabel, fontsize=20)
plt.title('Velocities relative to coma cluster vs projected distance',\
          fontsize=20)
plt.show()

#array of ecoids that fall within the shaded region (30 most likely real coma 
#candidates)
proper_eco = []
for index,value in enumerate(projected_distance):
    if value <= 2:
        if -1000 <= deltav.values[index] <= 1000:
            label = 'ECO'+labels[index]
            proper_eco.append(label)
          


          

            
            

    


