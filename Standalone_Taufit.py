# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa


#!/Applications/Spyder.app/Contents/MacOS/python


import argparse
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
args = parser.parse_args()


filepath = args.filename
pulseperiod = args.period

pulsar, nch, nbins, lm_rms = dri.read_header(filepath)

print "Pulsar name: %s" %pulsar
print "Number of channels: %d" %nch
print "Number of bins: %d" %nbins
print "RMS: %f" %lm_rms


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
    plt.plot(data,'y')
    plt.plot(noiselessmodel,'r')
    plt.title(r'$\tau: %.2f \pm %.4f$ sec' %(besttau*pulseperiod/nbins, taustd*pulseperiod/nbins))




# Fit power law through tau spectrum

# 1. Find error bars for tau

count = 200
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
        
        
#    x,y = np.histogram(e_besttaus)
#    plt.figure()
#    plt.hist(e_besttaus)
#    histobins = int(np.sqrt(count))
#    bincounts,binname = np.histogram(e_besttaus,bins=histobins)  ##Note tau here is in bins
#
#    sigmaline = int(25./100.*count) #quartiles
#    sum = 0
#    counter = 0
#    while (sum < sigmaline):
#        sum += bincounts[counter]
#        counter += 1
#    q1 = counter
#    sum = 0
#    counter = histobins - 1
#    while (sum < sigmaline):
#        sum += bincounts[counter]
#        counter -= 1
#    q3 = counter
#    
#    tau_low = binname[q1]
#    tau_high = binname[q3]
#
#    q1s.append(q1)
#    q3s.append(q3)
#    
#    print 'Histrogram Counter: %d' %counter

modelpow = Model(psr.PLM)

lmfittausstds = np.array(lmfittausstds)
obtainedtaus = np.array(obtainedtaus)
obtainedtausec = obtainedtaus*pulseperiod/nbins
freqms = np.array(freqmsMHz)/1000.

resultpow = modelpow.fit(e_besttaussec,x=e_freqms,weights=1,K=0.001,k=4)        
specfit = resultpow.best_values['k']
print 'alpha_1 = %.4f' %specfit

plt.figure()
plt.plot(1000*e_freqms,e_besttaussec,'ro')
plt.plot(1000*freqms,obtainedtausec,'c*',markersize=9.0)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)

plt.show()

