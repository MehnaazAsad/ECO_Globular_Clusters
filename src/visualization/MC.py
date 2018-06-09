#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: asadm2
"""
###DESCRIPTION
#This script runs an MC to get an estimate of number of GCs we could 
#potentially detect as well as an estimate of the halo mass given total mass
#of GC system. Histograms of 1000 realizations are plotted at the end.

import matplotlib.pyplot as plt
from scipy import integrate
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

#g-r colour of objects
g_r = []
for name in ECO_original.ECOID.values:
    gmag = ECO_photometric.gmag.loc[ECO_photometric.ECOID==name].values[0]
    rmag = ECO_photometric.rmag.loc[ECO_photometric.ECOID==name].values[0]
    colour = gmag-rmag
    g_r.append(colour)
mean_gr = np.mean(g_r)

###MC 
Mv_sun = 4.83 #V band magnitude of Sun
M_L_ratio = 1.45 #solar mass per solar luminosity
L_sun = 3.9*(10**33) #ergs/s
M_sun = 1.99*(10**33) #g
d = 97.6*(10**6) #pc

#Gaussian for luminosity distribution
def make_P_L(L_mu,sigma,L):
    constant = 1/(np.sqrt(2*np.pi)*sigma)
    e_num = (-2.5*np.log10(L/L_mu))**2
    e_denom = 2*(sigma**2)
    expression = constant*np.exp(-e_num/e_denom)
    return expression

#testing for ECO12028 currently but should eventually be run for all 
#single_halos
for name in single_halos.ECOID.values[:1]: 
    #g-band mag from ECO catalog
    gmag = ECO_photometric.gmag.loc[ECO_photometric.ECOID==name].values[0]
    #convert to M-band mag
    m_B  = gmag + (0.3130*(mean_gr)) + 0.2271 #Lupton transformation
    M_B = m_B-(5*np.log10(d - 5)) #apparent to absolute
    #calculate std in magnitude per galaxy
    sigma_g = 1.14 - (0.1*(M_B+20)) #from Virgo paper
    
    #parameters of gaussian for magnitude distribution
    mean_mag = -7.4 #from GC luminosity function in Virgo paper
    sigma_mag = sigma_g
    min_mag = -7.27 #calculated for this particlar galaxy based on distance, 
                    #redshift and ETC online
    
    #parameters of gaussian for number distribution
    mean_num = 10**(GC_count.log_GC_count.loc[GC_count.ECOID==name].values[0])
    sigma_num  = 0.3
    
    n_detected = []
    halo_masses = []
    #1000 realizations of MC
    for i in range(0,1000):
        #get number of GCs from gaussian
        num = np.random.normal(mean_num,sigma_num) 
        #assign magnitudes to num from gaussian
        mags = np.random.normal(mean_mag,sigma_mag,round(num)) 
        final_mags = [value for value in mags if value < min_mag] #mag cut
        n_detected.append(len(final_mags))
        
        #Convert magnitudes to luminosities
        gc_lum = []
        for value in final_mags:
            lum = 10**((value-Mv_sun)/(-2.5)) #per Lv_sun
            gc_lum.append(lum)
            
        mean_lum = 10**((mean_mag-Mv_sun)/(-2.5)) #per Lv_sun
        sigma_lum = sigma_mag            
        L_min = 10**((min_mag-Mv_sun)/(-2.5)) #per Lv_sun
        
        #Correction factor to account for rest of the GCs that are dimmer
        #than mag cut (turns out it doesn't make much of a difference since
        #bright GCs dominate light)
        L_obs = integrate.quad(lambda L: L*make_P_L(mean_lum,sigma_lum,L)\
                               ,L_min,10**10)
        L_tot = integrate.quad(lambda L: L*make_P_L(mean_lum,sigma_lum,L)\
                               ,0,10**10)
                             
        corr_factor = (L_obs[0]/L_tot[0])
        total_lum = sum(gc_lum)
        #Corrected observed total luminosity of GC system
        L_obs_corr = total_lum/corr_factor 
        #Total mass of GC system using M/L ratio
        M_GCS = np.log10(M_L_ratio*L_obs_corr)
        #USING GC_MASS-HALO_MASS RELATION
        M_h = np.log10(10**(M_GCS+4.15)) #M_sun
        halo_masses.append(M_h)
###
        
    #PLOTS
    #Figure 1: Histogram of number of GCs expected to detect (ECO12028)
    fig1 = plt.figure()
    n,bins,patches = plt.hist(n_detected,histtype = 'step')
    #Average of number estimate from MC
    x_peak = np.mean(n_detected)
    plt.axvline(x=x_peak,linestyle='dashed',color='b')
    plt.text(min(bins),max(n)-20,name,fontsize=15,\
             bbox=dict(facecolor='blue', alpha=0.5))
    plt.text(max(bins)-5,max(n)-20,x_peak,fontsize=15,\
             bbox=dict(facecolor='blue', alpha=0.5))
    plt.xlabel(r'Number of GCs')
    
    #Figure 2: Histogram of halo mass estimated for this object (ECO12028)
    fig2 = plt.figure()
    n,bins,patches = plt.hist(halo_masses,histtype = 'step')
    #Average of mass estimates from MC
    x_peak = float("%.2f" % round(np.mean(halo_masses),2))
    #Actual halo mass from ECO catalog
    x_ECO = single_halos.log_halogroupmass.loc[single_halos.ECOID==name]\
                .values[0]
    plt.axvline(x=x_peak,linestyle='dashed',color='b')
    plt.axvline(x=x_ECO,linestyle='dashed',color='y')
    plt.text(min(bins),max(n)-15,name,fontsize=15,\
             bbox=dict(facecolor='blue', alpha=0.5))
    plt.text(max(bins)-0.45,max(n)-15,x_peak,fontsize=15,\
             bbox=dict(facecolor='blue', alpha=0.5))
    plt.text(max(bins)-0.25,max(n)-15,x_ECO,fontsize=15,\
             bbox=dict(facecolor='yellow', alpha=0.5))
    plt.xlabel(r'$\log\ M_h/M_\odot$')