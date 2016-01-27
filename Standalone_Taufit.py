#!/usr/bin/env python

# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa


import argparse
import os, sys
import pypsr_redo as psr
import matplotlib.pyplot as plt
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
import numpy as np
from scipy import special
import DataReadIn as dri 

#"""Read and print header information"""

parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data files")
parser.add_argument('-p','--period',type=float,
                    help="Provide the period of the pulsar in seconds")
parser.add_argument('-c','--counts',type=int,
                    help="Provide the number of noise realisations (interger)")
args = parser.parse_args()


filepath = args.filename
pulseperiod = args.period
count = args.counts

newpath = r'./Plots' 
if not os.path.exists(newpath):
    os.makedirs(newpath)


pulsar, nch, nbins, lm_rms = dri.read_header(filepath)

print "Pulsar name: %s" %pulsar
print "Number of channels: %d" %nch
print "Number of bins: %d" %nbins
print "RMS: %f" %lm_rms

profilexaxis = np.linspace(0,pulseperiod,nbins)

#Fit profile using tau_fitter

obtainedtaus = []
lmfittausstds = []
freqmsMHz =[]
noiselessmodels =[]
comp_rmss = []

for i in range(nch):
    data, freqm = dri.read_data(filepath,i,nbins) 
    comp_rms = psr.find_rms(data,nbins)     
    noiselessmodel, besttau, taustd = psr.tau_fitter(data,nbins)
    
    obtainedtaus.append(besttau)    
    lmfittausstds.append(taustd)
    freqmsMHz.append(freqm)
    noiselessmodels.append(noiselessmodel)
    comp_rmss.append(comp_rms)

    
    plt.figure()
    plt.plot(profilexaxis,data,'y')
    plt.plot(profilexaxis,noiselessmodel,'r')
    plt.title('%s at %.1f MHz' %(pulsar, freqm))
    plt.text(0.55*pulseperiod,0.95*np.max(data),r'$\tau: %.4f \pm %.4f$ sec (lmfit error)' %(besttau*pulseperiod/nbins, taustd*pulseperiod/nbins),fontsize=14)    
    plt.ylim(ymax=1.1*np.max(data))    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('sec')
    
    plotname = '%s_Profile_fit_%d_%.1fMHz.png'  % (pulsar,i,freqm)
    print 'Saved %s in Plots' %plotname
    picpath = "/Users/marisa/Documents/PhD/GitHub_Scattering/BusyWeek/Plots"
    fileoutput = os.path.join(picpath,plotname)
    plt.savefig(fileoutput, dpi=150)



# Fit power law through tau spectrum

# 1. Find error bars for tau


#q1s = []
#q3s = []
#tau_lows = []
#tau_high = []

arr_ind = 0


e_freqms = np.zeros(count*nch)
e_besttaussec = np.zeros(count*nch)

for i in range(nch):
    rms = comp_rmss[i]
    profile =  noiselessmodels[i]
     
    for ii in range(count):
        noise = np.random.normal(0,rms,nbins)
        profile_noise = profile + noise
        e_noiselessmodel, e_besttaus, e_taustd = psr.tau_fitter(profile_noise,nbins)
        e_freqms[arr_ind] = freqmsMHz[i]/1000.
        e_besttaussec[arr_ind] = e_besttaus*pulseperiod/nbins
        arr_ind += 1
        print '%d/%d'% (arr_ind,count*nch)
        
        
##    x,y = np.histogram(e_besttaus)
##    plt.figure()
##    plt.hist(e_besttaus)
##    histobins = int(np.sqrt(count))
##    bincounts,binname = np.histogram(e_besttaus,bins=histobins)  ##Note tau here is in bins
##
##    sigmaline = int(25./100.*count) #quartiles
##    sum = 0
##    counter = 0
##    while (sum < sigmaline):
##        sum += bincounts[counter]
##        counter += 1
##    q1 = counter
##    sum = 0
##    counter = histobins - 1
##    while (sum < sigmaline):
##        sum += bincounts[counter]
##        counter -= 1
##    q3 = counter
##    
##    tau_low = binname[q1]
##    tau_high = binname[q3]
##
##    q1s.append(q1)
##    q3s.append(q3)
##    
##    print 'Histrogram Counter: %d' %counter
#
modelpow = Model(psr.PLM)

lmfittausstds = np.array(lmfittausstds)
obtainedtaus = np.array(obtainedtaus)
obtainedtausec = obtainedtaus*pulseperiod/nbins
freqms = np.array(freqmsMHz)/1000.

resultpowdata = modelpow.fit(obtainedtausec,x=freqms,weights=1,K=0.001,k=4)
resultpow = modelpow.fit(e_besttaussec,x=e_freqms,weights=1,K=0.001,k=4)    
specfitdata = resultpowdata.best_values['k']    
specfitspread = resultpow.best_values['k']

specdata_err= resultpowdata.params['k'].stderr
specspread_err = resultpow.params['k'].stderr

print 'from data: alpha = %.4f' %specfitdata
print 'from noise rms spread: alpha = %.4f' %specfitspread

plt.figure()
plt.plot(1000*e_freqms,e_besttaussec,'ro',label=r'$\tau$ spread from noise rms')
plt.plot(1000*freqms,obtainedtausec,'c*',markersize=9.0,label=r'$\tau$ from data')
plt.title(r'DATA: $\alpha = %.2f \pm %.2f$,  NOISE SPREAD:  $\alpha = %.2f \pm %.2f$' %(specfitdata,specdata_err,specfitspread,specspread_err), fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)
plt.legend()

plotnametau = 'Tau_fit_%s.png'  % (pulsar)
print 'Saved %s in Plots' %plotnametau
picpathtau = "/Users/marisa/Documents/PhD/GitHub_Scattering/BusyWeek/Plots"
fileoutputtau = os.path.join(picpathtau,plotnametau)
plt.savefig(fileoutputtau, dpi=150)


plt.show()

