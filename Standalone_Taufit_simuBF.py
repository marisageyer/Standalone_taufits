#!/usr/bin/python



# -*- coding: utf-8 -*-

#Created on Mon Jan 25 13:22:33 2016

#@author: marisa
#Standalone_Taufit_simu.py

import argparse
import os, sys
import pypsr_standalone as psr
import matplotlib.pyplot as plt
import lmfit
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
import numpy as np
from scipy import special
from scipy import stats
import DataReadIn as dri

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#"""Read and print header information"""
"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-f','--filename',
                    help="Provide the pathway to the data files")
parser.add_argument('-p','--period',type=float,
                    help="Provide the period of the pulsar in seconds")
parser.add_argument('-c','--counts',type=int,
                    help="Provide the number of noise realisations (interger)")
parser.add_argument('-s','--simulate',
                    help="Choosing this option leads to simulated data. For broadening function: choose between 'onedim', 'iso', 'aniso'.")                    
parser.add_argument('-m','--method',
                    help="Choosing method to fit data or simulation. Choose between 'onedim', 'iso', 'aniso','postfold', 'isoplusonedim'")                   
parser.add_argument('-r','--rawfile',
                    help="filepath to txt file, produced from shifting observing data")                   
parser.add_argument('-dc','--datacycle',
                    help="The type of data. Choose from comm., census, cycle5. Only used in creating filenames.")
parser.add_argument('-t','--template',
                    help="filepath to txt file, containing a high frequency EPN profile. It will use the fixed with of this profile.")



args = parser.parse_args()

"""Allocate variable names to the parsed options"""
filepath = args.filename
pulseperiod = args.period
count = args.counts
simu = args.simulate
meth = args.method
raw = args.rawfile
datac = args.datacycle
temp = args.template

"""Create folder to save to"""
newpath = r'./SummaryPlots_IsoOneDim'
if not os.path.exists(newpath):
    os.makedirs(newpath)


if simu is None:
    if raw is None:        
        print "\n Reading in data \n"
        pulsars = []
        pulsar, nch, nbins, lm_rms = dri.read_header(filepath)
        pulsars = [pulsar for i in range(nch)]
        print3 = "RMS: %f" %lm_rms
    elif raw.endswith(".ort.prof") == True:
        pulsar, nch, nbins, rawdata, freq = dri.read_Krishnak(raw)
        pulsars = [pulsar for i in range(nch)]
        print3 = ""
    elif raw.startswith('1937') == True or raw.startswith('../python-workingdir/Paul') == True:
        pulsar, nch, nbins, rawdata, freq = dri.read_Paul(raw)
        pulsars = [pulsar for i in range(nch)]
        print3 = ""
    else:
        pulsar, nch, nbins, rawdata, freq = dri.read_raw(raw)
        pulsars = [pulsar for i in range(nch)]
        print3 = ""
        
    print0 = "Pulsar name: %s" %pulsar
    print1 = "Number of channels: %d" %nch
    print2 = "Number of bins: %d" %nbins
    for k in range(4):
            print eval('print{0}'.format(k))
    
else:
    print "\n Simulating data \n"
    print(" \n 1. The %s broadening function" %simu)
    if simu in ('iso','ISO', 'Iso', 'onedim','1D','Onedim', 'postfold', 'pf', 'Postfold', 'POSTFOLD'):
        while True:
            try:
                taudivP = raw_input("Express max. tau as a fraction of the pulseperiod: ")
                taudivP = float(taudivP)
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break
    else:
        while True:
            try:
                taudivP1 = raw_input("Express max. tau1 as a fraction of the pulseperiod: ")
                taudivP1 = float(taudivP1)
                taudivP2 = raw_input("Express max. tau2 as a fraction of the pulseperiod: ")
                taudivP2 = float(taudivP2)
                taudivP = np.array([taudivP1, taudivP2])
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break
    print(" \n 2. The simulated pulsar")
    if pulseperiod is None:  
        while True:
            try:
                pulseperiod = raw_input("Express pulse period in seconds: ")
                pulseperiod = float(pulseperiod)
            except ValueError:
                print("Try again. Enter an int or float.")
                continue
            else:
                break  
    while True:
        try:
            dutycycle = raw_input("Express duty cycle as a percentage of pulse period, %s sec: " %pulseperiod)
            dutycycle = float(dutycycle)
        except ValueError:
            print("Try again. Enter an int or float.")
            continue
        else:
            break
    print(" \n 3. Data and observing properties")
    while True:
        try:
            freqlow, freqhigh, incr = raw_input("Choose a lowest and highest frequency and increment (MHz), e.g 50 120 5: ").split()
            freqlow, freqhigh, incr = float(freqlow), float(freqhigh), float(incr)
        except ValueError:
            print("Try again. Separate 3 floats with a space.")
            continue
        else:
            break

    while True:
        try:            
#            nbins = raw_input("Choose number of bins per pulse period: ")
#            nbins = int(nbins)
            snr = raw_input("Choose peak SNR: ")
            snr = int(snr)
        except ValueError:
            print("Try again. Enter an int or float.")
            continue
        else:
            break
    
    nbins = 512
    freqsimu = np.arange(freqlow,freqhigh+incr,incr)
    nch = len(freqsimu)
    propconst = taudivP*pulseperiod*(freqlow/1000.)**4
    propconst = np.array([propconst])
    
    tausecs = []
    for i in range(len(propconst)):
        tausec = propconst[i]*(freqsimu/1000.)**(-4)
        tausecs.append(tausec)
    tausecs = np.array(tausecs).transpose()
    
    pulsars = []
    for i in range(nch):
        if simu in ('aniso','Aniso','ANISO'):
            pulsar = r'Simul.: $\tau_1 = %.2f$, $\tau_2 = %.2f$ ' %(tausecs[i][0],tausecs[i][1])
            pulsars.append(pulsar)
        else:
            pulsar = r'Simul: $\tau_1 = %.2f$' %tausecs[i]
            pulsars.append(pulsar)
    pulsar = 'Simulated'
    print0 = "Pulsar : %s" %pulsar
    print1 = "Number of channels: %d" %nch
    print2 = "Number of bins: %d" %nbins
    print3 = ""

#Find pulseperiod from list if it wasn't parsed

pulsarBnamelist = ['B0037+56','B0114+58','B0540+23','B0611+22','B0740-28','B1848+12','B1907+10','B1911-04','B1915+13','B1920+21','B1933+16','B2255+58','B2303+30']
pulsarnamelist = ['J0040+5716','J0117+5914','J0543+2329','J0614+2229', 'J0742-2822','J1851+1259','J1909+1102','J1913-0440','J1917+1353','J1922+2110','J1935+1616','J2257+5909','J2305+3100']
pulsarperiodlist = [1.118225, 0.101439, 0.245975, 0.33496, 0.166762, 1.205303, 0.283641, 0.825936, 0.194631, 1.077924, 0.358738, 0.368246, 1.575886]


if pulseperiod == None: 
   if raw is not None and raw.endswith(".ort.prof") == True:
   	pulsarN = pulsar[0:8]
   	print pulsarN
   	pulseind = pulsarBnamelist.index(pulsarN)
   	pulseperiod = pulsarperiodlist[pulseind]
   else:
   	pulseind = pulsarnamelist.index(pulsar)
   	pulseperiod = pulsarperiodlist[pulseind]
else:
   pulseperiod = pulseperiod


profilexaxis = np.linspace(0,pulseperiod,nbins)

#for i in range(15):
#    plt.close(i)


## Create subplot structure
sp = 12 #number of subplots per figure
#plt.figure(figsize=(16,10))  

#Fit profile using tau_fitter

obtainedtaus = []
lmfittausstds = []
obtainedtaus2 = []
lmfittausstds2 = []

bestparamsall = []
bestparams_stdall = []
correls = []
redchis = []

freqmsMHz =[]
freqcsMHz =[]
noiselessmodels =[]
results = []
comp_rmss = []
comp_fluxes= []
comp_SNRs =[]
datas = []
climbvals =[]

besttau2 = []
taustd2 = []


#trainlength = 4

halfway = nbins/2.

for i in range(nch):
    if simu is None:
        if raw is None:
            data, freqc, freqm = dri.read_data(filepath,i,nbins)
            freqmsMHz.append(freqm)
            freqcsMHz.append(freqc)
            if i == 0:
                peakbin = np.argmax(data)
                print peakbin
            elif i==1 and pulsar in ('J1909+1102','J1913-0440') and datac in 'comm':
                peakbin = np.argmax(data)
                print peakbin
                print "I'm using the second channel to centre the data since the first one is dropped! This is for B1911-04 comm. data and B1907"
            else:
                peakbin = peakbin
                print peakbin
                shift = int(halfway-int(peakbin))
                print shift
            data = np.roll(data,int(halfway-int(peakbin)))
            print "I'm rolling the data to ensure the lowest freq. peak is in the middle"
        elif raw.endswith(".ort.prof") == True:
            data = rawdata
            freqmsMHz = freq*1000
        elif raw.startswith('1937') == True or raw.startswith('../python-workingdir/Paul') == True:
            data = rawdata
            freqmsMHz = freq*1000
        else:
            data = rawdata[i]
            freqmsMHz = freq*1000
    else:
        freqmsMHz = freqsimu
        freqGHz = freqmsMHz/1000.
        data = psr.simulate(pulseperiod,tausecs[i],dutycycle,-1.6,freqGHz[i],freqlow/1000.,nbins,snr,simu)
    comp_rms = psr.find_rms(data,nbins)
   
    if temp is None: 
        if meth is None:
            print "No fitting method was chosen. Will default to an isotropic fitting model. \n Use option -m with 'onedim' or 'aniso' to change."
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_fitter(data,nbins)
            #bestparams_order: sig, mu, A, DC
            climbval = psr.returnclimb(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        elif meth in ('iso','Iso','ISO'):
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_fitter(data,nbins)
            climbval = psr.returnclimb(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        elif meth in ('onedim','1D','Onedim'):
            result, noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi, corsig = psr.tau_1D_fitter(data,nbins)
            climbval = psr.returnclimb1D(np.linspace(1,nbins,nbins),bestparams[1],bestparams[0],bestparams[2],besttau,bestparams[3],nbins)
        elif meth in ('aniso','Aniso','ANISO'):
            noiselessmodel, besttau, taustd, besttau2, taustd2, bestparams, bestparams_std, redchi = psr.tau_ani_fitter(data,nbins)
        elif meth in ('fix'):
              #This uses fixed width values, as given by this fixval array
              # I have also used it to do fixed taus - double check the current state of psr.tau_fitter_fix before using it
#              fixvals = np.array([ 166.79287704,   96.42112647,   43.33346859,   31.37441266,
#         22.98447972,   20.0436986 ,   21.36908006,   22.84946732,
#         19.17358203,   18.54216172,   15.88792868,   14.93648958,
#         12.94439522,   11.45639279,   13.1817301 ,   10.61473998])
              fixvals = np.array([ 735.96503772,  570.80448298,  283.9382534 ,  220.03587608,
        157.5615805 ,  142.65554165,  106.30809045,  114.00090536])
#              fixvals = np.array([4.03238179,  3.7952635 ,  3.57981438,  3.38564806,  3.21046153, 3.04428873,  2.89688641,  2.76022124,  2.63403959,  2.51735798, 2.40919727,  2.30870483,  2.21432137, 2.12914729,  2.0447208 ,1.97032497])
              noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fitter_fix1D(data,nbins,fixvals[i])
        elif meth in ('multi'):
              if i == 0:
#        		    bw1 = raw_input("Provide bin-window for 1st the peak (e.g: [20,30]): ")
#        		    bw1 = eval(bw1)
#        		    bw2 = raw_input("Provide bin-window start and end 2nd peak (e.g: [35,45]): ")
#        		    bw2 = eval(bw2)
        		    bw1, bw2 = np.linspace(660,720,2), np.linspace(725,750,2)        
        		    peak1 = np.max(data[bw1[0]: bw1[-1]])
        		    peak2 = np.max(data[bw2[0]: bw2[-1]])
        		    maxmean = 2.0*bw2[-1]
              else:
                     bw1,bw2 = bw1,bw2
                     peak1, peak2, maxmean = peak1, peak2, maxmean
              nms, bt, tstd, bp, bpstd, chis = [], [], [], [], [], []
              for k in range(0,len(bw1)):
                  for j in range(0,len(bw2)):
        			noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_fittermulti(data,nbins,bw1[k],bw2[j], peak1, peak2, maxmean)
        			nms.append(noiselessmodel)
        			bt.append(besttau)
        			tstd.append(taustd)
        			bp.append(bestparams)
        			bpstd.append(bestparams_std)
        			chis.append(redchi)
        			print j*k
              tstd = np.array(tstd)
              nonz = np.nonzero(tstd)
              nonz_min = np.argmin(tstd[np.nonzero(tstd)])
              tstd_ind = np.array(nonz).flatten()[nonz_min]     
              noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = nms[tstd_ind], bt[tstd_ind], tstd[tstd_ind], bp[tstd_ind], bpstd[tstd_ind], chis[tstd_ind]               
        elif meth in ('postfold', 'pf', 'Postfold', 'POSTFOLD'):
              trainlength = np.ceil(8*tausecs/pulseperiod)
              trainlength = np.array([int(x) for x in trainlength])
              print trainlength
              noiselessmodel, besttau, taustd = psr.tau_fitter_postfold(data,nbins,trainlength[i])
        else:
             print "Incorrect fitting method. Choose from iso, onedim, aniso, postfold"
    else:
        templatefile = np.loadtxt(temp)
        if meth is None:
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_tempfitter(data,nbins,tempdata)
        elif meth in ('iso','Iso','ISO'):
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_tempfitter(data,nbins,tempdata)
        elif meth in ('onedim','1D','Onedim'):
            tempdata = templatefile[:,3]
            noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = psr.tau_1D_tempfitter(data,nbins,tempdata)                
        else:
           print "Incorrect fitting method. Choose from iso or onedim or multi (for analytical double peak)"
    


    """ INPUT SECTION ENDS """

            
    comp_flux = psr.find_modelflux(noiselessmodel,nbins)
    comp_SNR_model = psr.find_peaksnr(noiselessmodel,comp_rms)

    print 'SNR (from model): %.2f' % comp_SNR_model
    comp_SNR =  psr.find_peaksnr_smooth(data,comp_rms)
    print 'SNR (from data): %.2f' % comp_SNR
#    print "The SNR is computed using data only, not model"
#    print "The SNR is computed using model"
    
    obtainedtaus.append(besttau)    
    lmfittausstds.append(taustd)
    obtainedtaus2.append(besttau2)    
    lmfittausstds2.append(taustd2)
    bestparamsall.append(bestparams)
    bestparams_stdall.append(bestparams_std)
    redchis.append(redchi)
    correls.append(corsig)
    
        
    noiselessmodels.append(noiselessmodel)
    results.append(result)
    comp_rmss.append(comp_rms)
    comp_fluxes.append(comp_flux)
    comp_SNRs.append(comp_SNR)
    datas.append(data)
    climbvals.append(climbval)


"Pick out the correlations of sigma with Amp to use in flux errors"

cor_sigA = np.zeros(len(correls))
for i in range(len(correls)):
    if correls[i] is None:
        print i
        cor_sigA[i] = 0
    elif correls[i] is not None:
        cor_sigA[i] = correls[i]['A']
#        cor_sigA.append(csigA)

"""Considering the generated BW, calculate the monochromatic frequency that should be associated with the tau-fits"""
"""Replace freqmsMHz with this frequency"""
#for now I'm only using an average bandwidth for all observations. Considering changing it from freq to freq.

"""Fix SNR cutoffs, include couple of specific ones for individual pulsars"""
	
SNRcutoff = 2.65
tauerrorcutoff = 500.0

if pulsar in 'J1913-0440' and datac in 'comm':
    SNRcutoff = 6.5
    

if nch !=1:
    """Delete profiles with 

    1. inadequate SNR values,
    2. too large tau error values
    3. tau error == 0, i.e. not a fit

    update other arrays accordingly"""

    print4 = "SNR cutoff: %.2f" % SNRcutoff
    print5 = "Tau perc error cutoff: %.2f" % (100*tauerrorcutoff)

    for k in range(0,6):
        print eval('print{0}'.format(k))
            
    if meth in 'fix':
        ##This is for when tau is fixed, because the errors in lmfittau is then 0 causing the code to crash
        data_highsnr = np.array(datas)
        model_highsnr = np.array(noiselessmodels)
        taus_highsnr = np.array(obtainedtaus)
        lmfitstds_highsnr = np.array(lmfittausstds)
        taus2_highsnr = np.array(obtainedtaus2)
        lmfitstds2_highsnr = np.array(lmfittausstds2)
        freqMHz_highsnr = np.array(freqmsMHz)
        freqms_highsnr = np.array(freqMHz_highsnr)/1000.
        fluxes_highsnr = np.array(comp_fluxes)
        rms_highsnr = np.array(comp_rmss)
        ind_lowSNR = (np.array([],),)
        ind_tauerr = (np.array([],),)

        
    else:        
        ind_lowSNR = np.where(np.array(comp_SNRs) < SNRcutoff)
        ind_tauerr = np.where(np.array(lmfittausstds)/np.array(obtainedtaus) > tauerrorcutoff)
        ind_tauerr2 = np.where(np.array(lmfittausstds)==0)
        ind_tauERRs = np.hstack((ind_tauerr,ind_tauerr2))
        ind_uniq = np.unique(np.hstack((ind_lowSNR,ind_tauerr,ind_tauerr2)))
        data_highsnr = np.delete(np.array(datas),ind_uniq,0)
        model_highsnr = np.delete(np.array(noiselessmodels),ind_uniq,0)
    
        """"""
    
        taus_highsnr = np.delete(np.array(obtainedtaus),ind_uniq)
        lmfitstds_highsnr = np.delete(np.array(lmfittausstds),ind_uniq)
        taus2_highsnr = np.delete(np.array(obtainedtaus2),ind_uniq)
        lmfitstds2_highsnr = np.delete(np.array(lmfittausstds2),ind_uniq)
        climb_highsnr = np.delete(np.array(climbvals),ind_uniq)
    
    
        freqMHz_highsnr = np.delete(np.array(freqmsMHz),ind_uniq)
        freqms_highsnr = np.array(freqMHz_highsnr)/1000.
        fluxes_highsnr = np.delete(np.array(comp_fluxes),ind_uniq)
        corsigA_highsnr = np.delete(np.array(cor_sigA),ind_uniq)
        rms_highsnr = np.delete(np.array(comp_rmss),ind_uniq)
        

    number_of_plotted_channels = len(data_highsnr)
    npch = number_of_plotted_channels

    """Array with all the other fitting parameters: sigma, A, etc."""
    bestpT = np.transpose(bestparamsall)
    bestpT_std = np.transpose(bestparams_stdall)

    print6 = "Number of plotted channels: %d/%d" %(npch, nch)
    print7 = "Number of channels dropped for low SNR: %d" %np.shape(ind_lowSNR)[1]
    print8 = "Number of channels dropped for high tau error or zero tau error: %d" %np.shape(ind_tauERRs)[1]

    for k in range(6,9):
        print eval('print{0}'.format(k))

    """Similarly reduce array to only the used data (i.e. lowSNR and hightauerr removed)"""
    bestpT_highSNR = np.zeros([len(bestpT),npch])
    bestpT_std_highSNR = np.zeros([len(bestpT),npch])
    if meth in 'fix':
        bestpT_highSNR = bestpT
        bestpT_std_highSNR = bestpT_std
    else:
        for i in range(len(bestpT)):
            bestpT_highSNR[i]= np.delete(bestpT[i],ind_uniq)
            bestpT_std_highSNR[i]= np.delete(bestpT_std[i],ind_uniq)

    """Calculate fits for parameters sigma and mu"""
    """Shift data based on mu-trend (DM trend?)"""

    pbs = pulseperiod/nbins
    


    """Fit models to sigma"""
    powmod = PowerLawModel()
    powpars = powmod.guess(bestpT_highSNR[0], x=freqms_highsnr)
    powout = powmod.fit(bestpT_highSNR[0], powpars, x=freqms_highsnr, weights=1/((bestpT_std_highSNR[0])**2))

    linmod = LinearModel()

    quadmod = QuadraticModel()
    quadpars = quadmod.guess(bestpT_highSNR[0], x=freqms_highsnr)
    quadout  = quadmod.fit(bestpT_highSNR[0], quadpars, x=freqms_highsnr, weights=1/((bestpT_std_highSNR[0])**2))

    expmod = ExponentialModel()
    exppars = expmod.guess(bestpT_highSNR[0], x=freqms_highsnr)
    expout = expmod.fit(bestpT_highSNR[0], exppars, x=freqms_highsnr, weights=1/((bestpT_std_highSNR[0])**2))

    """Fit quadractic/DM model to mu"""
#    quadpars_mu = quadmod.guess(bestpT_highSNR[1], x=freqms_highsnr)
#    quadout_mu  = quadmod.fit(bestpT_highSNR[1], quadpars_mu, x=freqms_highsnr, weights=1/((bestpT_std_highSNR[1])**2))

    """Fit a DM model to delta mu"""
    delnuarray = [-(1/freqMHz_highsnr[-1]**2-1/freqMHz_highsnr[i]**2) for i in range(npch)] ##in MHz
    delmuarray = [(bestpT_highSNR[1][-1] - bestpT_highSNR[1][i])*pbs for i in range(npch)] ##in seconds
    delmu_stdarray = [(bestpT_std_highSNR[1][-1] - bestpT_std_highSNR[1][i])*pbs for i in range(npch)]

    DM_linpars = linmod.guess(delmuarray, x=delnuarray)
#	DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray, weights=1/(np.power(delmu_stdarray,2)))
    DM_linout  = linmod.fit(delmuarray, DM_linpars, x=delnuarray)

    DM_CCval = DM_linout.best_values['slope']
    DM_CCvalstd = DM_linout.params['slope'].stderr

    DMmodelfit = DM_linout.best_fit ##model gives deltime in seconds (used to shift data)

    DMconstant = 4148.808
    #uncertainty in the constant is 0.003 - only affects the Delta DM value in the 9th decimal
    DMval = (DM_CCval/DMconstant)
    DMvalstd = (DM_CCvalstd/DMconstant)
    DMcheck = psr.DM_checker(freqmsMHz,bestpT_highSNR[1]*pbs)


    ###Calculate the required shift in bins as an integer
    mu_shift_int = np.around(DMmodelfit/pbs).astype(int)

    #if raw is None:
    #    shiftpath = r'./ShiftedTxtfiles'
    #    if not os.path.exists(shiftpath):
    #        os.makedirs(shiftpath)
    #    shiftfile = '%s_rawshiftdata_%s.txt' %(pulsar,meth)
    #    shiftpathfile = os.path.join(shiftpath,shiftfile)   
    #
    #    shifted_data = [np.roll(datas[i],mu_shift_int[i]) for i in range(npch)]
    #    freqsshape = freqms_highsnr.reshape(npch,1)
    #    shifted_stackfile = np.hstack((shifted_data,freqsshape))
    #    np.savetxt(shiftpathfile,shifted_stackfile)

else:
    plt.figure(1,figsize=(6,4))
    plt.rc('text', usetex = True)
    plt.rc('font',family='serif')
    plt.plot(data,'m',alpha=0.8, label="data")
    plt.plot(noiselessmodel, 'c', alpha=0.7, label="model")
    #unscat = psr.TwoPeaksModel(np.linspace(1,nbins,nbins), 1024, bp[tstd_ind][4], bp[tstd_ind][5], bp[tstd_ind][2] , bp[tstd_ind][3], bp[tstd_ind][0], bp[tstd_ind][1], bp[tstd_ind][6])
    #plt.plot(unscat,'g--')    
    plt.title('PSR %s at %.1f MHz' %(pulsars[0], freqmsMHz))
    plt.text(nbins/1.5,0.8*np.max(data),r'$\tau: %.4f \pm %.1e$ms' %(besttau*pulseperiod/(2*nbins)*1000, taustd*pulseperiod/(2*nbins)*1000),fontsize=11)
    plt.legend()    
#    plt.xlim(xmin=500,xmax=1024)
#    plt.annotate(r'$\tau: %.4f \pm %.4f$ ms' %(besttau*pulseperiod/nbins*1000, taustd*pulseperiod/nbins*1000),xy=(np.max(profilexaxis),np.max(data)),xycoords='data',xytext=(0.8,0.8),textcoords='axes fraction',fontsize=12)

#    Kplot = '%s_%s.png'  % (pulsar,raw[-6:-4])
#    picpathtau = newpath
#    fileoutputtau = os.path.join(picpathtau,Kplot)
#    plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s.txt' % (pulsar,raw[-6:-4]), np.column_stack((bt[tstd_ind], tstd[tstd_ind])).flatten())
#    plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s_params.txt' % (pulsar,raw[-6:-4]), bp[tstd_ind])
#    plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s_paramstd.txt' % (pulsar,raw[-6:-4]), bpstd[tstd_ind])
#    print 'Saved %s in %s' %(Kplot,newpath) 
    sys.exit()
    

    
"""Plotting starts"""

if meth in 'onedim':
    print "METH IS ONEDIM"
    alfval = 0.6
    alfval2 = 0.2
    markr = 'k^'
    markr2 = 'b^'
    prof = 'b--'
    ls = 'dashed'
    lcol='b'
    textpos = 0.7
    textpos2 = 3
if meth in 'iso':
    print "METH IS ISO"
#    plt.close('all')
    alfval = 1.0
    alfval2 = 0.4
    markr = 'k*'
    markr2 = 'r*'
    prof = 'r-'
    ls='solid'
    lcol ='r'
    textpos = 0.8
    textpos2 = 5
else:
    print "NO METH SPECIFIED"
    
##PLOT PROFILES##

totFig = (npch+5)/sp + 1
#for i in range(totFig):
#    plt.figure((totFig-i), figsize=(16,10))
#    plt.figure((totFig-i))

profiles = []
for j in range(npch):
    if j+1 == sp:
        numFig = (j+1)/sp
        subplotcount = sp
    else: 
        numFig = (j+1)/sp + 1
        if j+1 < sp:
            subplotcount = j+1
        else: 
            subplotcount = j+1 - sp
    bins, profile = psr.makeprofile(nbins = nbins, ncomps = 1, amps = bestpT_highSNR[2][j], means = bestpT_highSNR[1][j], sigmas = bestpT_highSNR[0][j])
    profiles.append(profile)
    smootheddata = psr.smooth(datas[j],int(0.05*nbins))
    figg = plt.figure(numFig,figsize=(16,10))
    figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
    plt.subplot(3,4,subplotcount)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
#    plt.tight_layout()
#    data_highsnr[j] = np.hstack((data_highsnr[j][800:1024],data_highsnr[j][0:800]))
#    model_highsnr[j] = np.hstack((model_highsnr[j][800:1024],model_highsnr[j][0:800]))
    plt.plot(profilexaxis,data_highsnr[j],'y',alpha = 0.25)
    plt.plot(profilexaxis,model_highsnr[j],prof, alpha = 0.7)
#    plt.plot(profilexaxis,shifted_data[j],'g',alpha = 0.6)
#    plt.plot(profilexaxis,smootheddata,'g', alpha = 0.4,lw=2.0)
#    plt.plot(profilexaxis,profiles[j],'k', alpha = 0.7)
    plt.title('PSR %s at %.1f MHz' %(pulsars[j], freqMHz_highsnr[j]))
    plt.annotate(r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[j]*pulseperiod/nbins, lmfitstds_highsnr[j]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[j])),xycoords='data',xytext=(0.4,textpos),textcoords='axes fraction',fontsize=12)
#    if len(taus2_highsnr) != 0:
#        plt.annotate(r'$\tau2: %.5f \pm %.5f$ sec' %(taus2_highsnr[j]*pulseperiod/nbins, lmfitstds2_highsnr[j]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[j])),xycoords='data',xytext=(0.4,0.6),textcoords='axes fraction',fontsize=12)        
    plt.ylim(ymax=1.1*np.max(data_highsnr[j]))
#    plt.ylim(ymax=1.1*np.max(profiles[j]))
#    plt.ylim(-200.0,700.0)
    plt.xlim(xmax=pulseperiod)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('time (sec)')
    plt.ylabel('intensity (mJy)')
        

if npch >= sp:
    subpl_ind = int(npch-sp) +1
else:
    subpl_ind = int(npch) + 1

if npch + 1 == sp:
    numfig = (npch+1)/sp
else:
    numfig = (npch+1)/sp + 1

##PLOT FLUX##

unscatflux = []
for i in range(npch):
    unscatfl = np.sum(profiles[i])/nbins
    unscatflux.append(unscatfl)

"""Calculate error in Flux"""

sigmaWIDTH = bestpT_std_highSNR[0]*pbs #in seconds
sigmaAMP = bestpT_std_highSNR[2]  #in mJy
WIDTHS =bestpT_highSNR[0]*pbs #in seconds
AMPS =bestpT_highSNR[2] #in mJy

Expr1 = np.sqrt(2*np.pi)*AMPS
Expr2 = np.sqrt(WIDTHS)
AreaExpression = Expr1*Expr2

sigmaEx1 = np.sqrt(2*np.pi)*sigmaAMP
sigmaEx2 = Expr2*0.5*sigmaWIDTH/WIDTHS

sigmaFlux =AreaExpression*np.sqrt(np.power(sigmaEx1/Expr1,2)+ np.power(sigmaEx2/Expr2,2)+2*corsigA_highsnr*sigmaEx1*sigmaEx2/(Expr1*Expr2))

correctedflux = fluxes_highsnr+climb_highsnr
meancorflux = np.mean(correctedflux)
meancorfluxerr = np.sqrt(np.sum(correctedflux**2))/len(correctedflux)

plt.figure(numfig, figsize=(16,10))
plt.subplot(3,4,subpl_ind)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(freqMHz_highsnr,fluxes_highsnr,alpha = alfval)
#plt.plot(freqMHz_highsnr, fluxes_highsnr+bestpT_highSNR[3],prof,markersize=9.0,linewidth=1.0)
plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=0.3)
plt.errorbar(freqMHz_highsnr, unscatflux,yerr=sigmaFlux, fmt=markr2,markersize=10.0, alpha=alfval2)
plt.title('PSR %s' %(pulsar))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
plt.ylabel(r'Calibrated flux (mJy)',fontsize=16)

lmfittausstds = np.array(lmfittausstds)
obtainedtaus = np.array(obtainedtaus)
lmfittausstds2 = np.array(lmfittausstds2)
obtainedtaus2 = np.array(obtainedtaus2)

obtainedtausec = obtainedtaus*pulseperiod/nbins
obtainedtausec2 = obtainedtaus2*pulseperiod/nbins
lmfitstdssec_highsnr = lmfitstds_highsnr*pulseperiod/nbins
lmfitstdssec2_highsnr = lmfitstds2_highsnr*pulseperiod/nbins

freqms = np.array(freqmsMHz)/1000. #(Use this one to compare the alpha outcome of all the freq channels with only the high SNR/low tau error ones.)

"""The associated high SNR arrays are"""
taus_highsnr = np.array(taus_highsnr)
taus2_highsnr = np.array(taus2_highsnr)
taussec_highsnr = taus_highsnr*pulseperiod/nbins
taussec2_highsnr = taus2_highsnr*pulseperiod/nbins


"""Order taus in the case of an anisotropic model with 2 taus"""
taus_order, taus2_order = np.zeros(len(taus_highsnr)), np.zeros(len(taus_highsnr))
taus_order_err, taus2_order_err = np.zeros(len(taus_highsnr)), np.zeros(len(taus_highsnr))

if len(taussec2_highsnr) != 0:
    for i in range(len(taus_highsnr)):
        if taussec_highsnr[i] > taussec2_highsnr[i]:
            taus_order[i] = taussec_highsnr[i]
            taus2_order[i] = taussec2_highsnr[i]
            taus_order_err[i] = lmfitstdssec_highsnr[i]
            taus2_order_err[i] = lmfitstdssec2_highsnr[i]     
        else:
            taus_order[i] = taussec2_highsnr[i]
            taus2_order[i] = taussec_highsnr[i]
            taus_order_err[i] = lmfitstdssec2_highsnr[i]
            taus2_order_err[i] =  lmfitstdssec_highsnr[i]
else: 
    taus_order = taussec_highsnr

if len(taussec2_highsnr) != 0:
    taussec_highsnr = taus_order
    taussec2_highsnr = taus2_order
    tausgeo = np.sqrt(taussec_highsnr *taussec2_highsnr)


"""CALCULATE FITS TO TAUS"""
powmod = PowerLawModel()
#powparstau = powmod.guess(obtainedtausec,x=freqms)
#powouttau = powmod.fit(obtainedtausec,powparstau, x=freqms,weights=1/(lmfittausstds**2))

powparstau_highsnr = powmod.guess(taussec_highsnr,x=freqms_highsnr)
powouttau_highsnr = powmod.fit(taussec_highsnr,powparstau_highsnr, x=freqms_highsnr,weights=1/(lmfitstdssec_highsnr**2))

#specfitdata = -powouttau.best_values['exponent']
specfitdata_highsnr = -powouttau_highsnr.best_values['exponent']
specdata_err_highsnr = powouttau_highsnr.params['exponent'].stderr

spec_amp = powouttau_highsnr.best_values['amplitude']
spec_err_amp = powouttau_highsnr.params['amplitude'].stderr

tau_HBAtop = spec_amp*np.power(freqms_highsnr[-1],-specfitdata_highsnr)
tau_1GHz = psr.tauatfreq(freqms_highsnr[-1],tau_HBAtop,1.0,specfitdata_highsnr)
tau_100MHz = psr.tauatfreq(freqms_highsnr[-1],tau_HBAtop,0.1,specfitdata_highsnr)

#print ""
#print9 = 'alpha (from all data) = %.4f' %specfitdata
print9 =""
print10 = 'alpha (from high SNR, low tau err) = %.4f' %specfitdata_highsnr
print11 = 'pulseperiod = %.6f' %pulseperiod

print ""
print freqcsMHz
print freqmsMHz
print ""
#print freqms_highsnr
#print freqMHz_highsnr

for k in range(9,12):
    print eval('print{0}'.format(k))

print "tau_1GHz = %.8f sec" %tau_1GHz
print "tau_100MHz = %.8f sec" %tau_100MHz
    
##PLOT TAU##  
    
if subpl_ind >= sp:
    subpl_ind2 = int(subpl_ind-sp) +1
else:
    subpl_ind2 = int(subpl_ind) + 1

if npch + 2 == sp:    
    numfig2 = (npch+2)/sp
else:
    numfig2 = (npch+2)/sp + 1

plt.figure(numfig2, figsize=(16,10))
plt.subplot(3,4,subpl_ind2)
#plt.tight_layout()
plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,markersize=9.0,capthick=2,linewidth=1.5, alpha = alfval)
plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-', alpha=alfval)
#plt.plot(1000.*freqms_highsnr, powouttau_highsnr1.best_fit, 'k--', label=r'Unweighted: $\alpha = %.2f \pm %.2f$' %(specfitdata_highsnr1,specdata_err_highsnr1))
plt.title(r'$\alpha = %.2f \pm %.2f$' %(specfitdata_highsnr,specdata_err_highsnr))
#plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)



### PLOT TAU LOG ##  

if subpl_ind2 >= sp:
    subpl_ind3 = int(subpl_ind2-sp) +1
else:
    subpl_ind3 = int(subpl_ind2) + 1
    
if npch + 3 == sp:    
    numfig3 = (npch+3)/sp
else:
    numfig3 = (npch+3)/sp + 1


freqmsMHz = np.array(freqmsMHz)
ticksMHz = freqmsMHz.astype(np.int)[0:len(freqmsMHz):2]

plt.figure(numfig3,figsize=(16,10))
plt.subplot(3,4,subpl_ind3)
#plt.figure(NF,figsize=(12,4))
#plt.subplot(1,3,1)
plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,alpha = alfval, markersize=11.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(specfitdata_highsnr,specdata_err_highsnr))  
plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-', alpha = alfval)
plt.title('%s' %pulsar)
#plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)
plt.legend(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim(xmin=(freqms[0])*950,xmax=(freqms[-1])*1050)
plt.xticks(ticksMHz,ticksMHz,fontsize=12)
#plt.tight_layout()


##PLOT CHANGE in ALPHA ##  

if subpl_ind3 >= sp:
    subpl_ind4 = int(subpl_ind3-sp) +1
else:
    subpl_ind4 = int(subpl_ind3) + 1
    
if npch + 4 == sp:    
    numfig4 = (npch+4)/sp
else:
    numfig4 = (npch+4)/sp + 1

"""Calculate how the spectral index is changing as more and more low freq. channels are added"""

spec_sec = np.zeros(npch-1)
spec_std_sec = np.zeros(npch-1)

for i in range(npch-1):
    bb = npch - (i+2)
    ee = npch
#    print freqms[bb:ee]
    powparstau_sec = powmod.guess(taussec_highsnr[bb:ee],x=freqms_highsnr[bb:ee])
    powouttau_sec = powmod.fit(taussec_highsnr[bb:ee],powparstau_sec, x=freqms_highsnr[bb:ee],weights=1/(lmfitstds_highsnr[bb:ee]**2))    
    spec_sec[i] = -powouttau_sec.best_values['exponent']
    spec_std_sec[i] = powouttau_sec.params['exponent'].stderr

freq_incl = 1000*freqms_highsnr[::-1][1:]

plt.figure(numfig4, figsize=(16,10))
plt.subplot(3,4,subpl_ind4)
plt.errorbar(freq_incl,spec_sec,yerr=spec_std_sec,fmt=markr,alpha=alfval,markersize=6.0,capthick=2,linewidth=0.5)
plt.title(r'Change in spectral index')
for x,y in zip(freq_incl[0::2], spec_sec[0::2]):
    plt.annotate('%.2f' %y, xy=(x,1.08*y), xycoords='data',textcoords='data')
plt.ylim(plt.ylim()[0],1.1*plt.ylim()[1])
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel(r'Lowest freq included (MHz)',fontsize=16)
plt.ylabel(r'$\alpha$',fontsize=16)
#plt.tight_layout()

#
###PLOT CHI ##  

if subpl_ind4 >= sp:
    subpl_ind5 = int(subpl_ind4-sp)+1
else:
    subpl_ind5 = int(subpl_ind4)+1

if npch + 5 == sp:    
    numfig5 = (npch+5)/sp
else:
    numfig5 = (npch+5)/sp + 1

redchis_highsnr = np.delete(np.array(redchis),ind_uniq)

plt.figure(numfig5, figsize=(16,10))
plt.subplot(3,4,subpl_ind5)
plt.plot(1000.*freqms_highsnr, redchis_highsnr/np.power(rms_highsnr,2), markr,alpha=alfval,markersize = 12)
#cres = []
#for i in range(npch):
#    varnr = 5
#    compres = (np.sum(np.power(data_highsnr[i] - model_highsnr[i],2))/(nbins-varnr))/np.power(rms_highsnr[i],2)
#    cres.append(compres)
#plt.plot(1000.*freqms_highsnr, cres, 'm*', markersize = 12)
plt.title(r'Reduced $\chi^2$ values')
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=10)
plt.xticks(fontsize=12)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
plt.ylabel(r'$\chi^2$',fontsize=16)

"""Compute the KS Test D value and probability of the residuals following a normal distribution"""
KSs = np.zeros((npch,2))
ADs = np.zeros((npch,2))
for i in range(npch):
    resdata = data_highsnr[i] - model_highsnr[i]
    resnormed = (resdata-resdata.mean())/resdata.std()
    KSd, KSp = stats.kstest(resnormed, 'norm')
    KSs[i,0] = KSd
    KSs[i,1]= KSp
    print KSd, KSp
    
#    aa,bb,cc = stats.anderson(resnormed, 'norm')
#    print aa,bb,cc

resmeanP = np.mean(KSs[:,1])
print "Mean probability of residuals being Gaussian: %.4f" % resmeanP

print "Mean reduced Chi square: %.2f" %  np.mean(redchis_highsnr/np.power(rms_highsnr,2))

for i in range(npch):
    print "Corrected Flux: %.1f MHz, %.1f \pm %.1f" %(freqms_highsnr[i]*1000,correctedflux[i],sigmaFlux[i])

print "Max flux:"
print np.max(correctedflux + sigmaFlux)
print "Min flux:"
print np.min(correctedflux - sigmaFlux)
print "Mean flux"
print np.mean(correctedflux)
print "DM"
print "%.4f pm %.4f" %(DMval,DMvalstd)

#### PLOT MU ##
#
#plt.figure(numfig5)
#plt.subplot(3,4,subpl_ind5)
##plt.errorbar(bestpT_highSNR[1]*pbs,freqMHz_highsnr, xerr = bestpT_std_highSNR[1]*pbs, fmt = 'b*',markersize=9.0,capthick=2,linewidth=1.5)
##plt.plot(quadout_mu.best_fit*pbs, freqMHz_highsnr, 'b-', label='a,b = %.2f,%.2f' %(quadout_mu.best_values['a'],quadout_mu.best_values['b']))
#plt.errorbar(np.array(delmuarray),1/freqMHz_highsnr, xerr = delmu_stdarray, fmt = 'b*',markersize=9.0,capthick=2,linewidth=1.5)
##plt.plot(DMbestfit, 1/freqMHz_highsnr, 'b-', label=r'DM = %.4f \pm %.6f' %(DMmyval,DMmyval_std))
##plt.plot(1000.*freqms_highsnr, powout_mu.best_fit*pbs, 'm-', label='pow = %.2f' %powout_mu.best_values['exponent'])
##plt.plot(quadout_mu.best_fit*pbs,1000.*freqms_highsnr, 'b-', label='a,b = %.2f,%.2f' %(quadout_mu.best_values['a'],quadout_mu.best_values['b']))
##plt.plot(bestpT_highSNR[1]*pbs,1000*freqms_highsnr, 'k*-', label=r'DM: %.4f pc.cm^{-3}' %(DMval))
#plt.yticks(fontsize=10)
#plt.xticks(fontsize=7)
#plt.xlabel(r'$\mu (sec)$ MHz',fontsize=16)
#plt.ylabel(r'$\nu$ MHz',fontsize=16)
#plt.legend(fontsize = 9)

#powparstau_highsnr = powmod.guess(taussec_highsnr,x=freqms_highsnr)
#powouttau_highsnr = powmod.fit(taussec_highsnr,powparstau_highsnr, x=freqms_highsnr,weights=1/(lmfitstdssec_highsnr**2))



if len(taussec2_highsnr) != 0:
    resultpow2guess = powmod.guess(taussec2_highsnr,x=freqms_highsnr)
    resultpowdata2_highsnr = powmod.fit(taussec2_highsnr,resultpow2guess, x=freqms_highsnr,weights=1/(lmfitstdssec2_highsnr**2),K=0.001,k=4)
    
    resultpowgeoguess = powmod.guess(tausgeo,x=freqms_highsnr)
    resultpowdatageo_highsnr = powmod.fit(tausgeo,resultpowgeoguess, x=freqms_highsnr,weights=1/(lmfitstdssec2_highsnr**2),K=0.001,k=4)
    
    specfitdata2_highsnr = -resultpowdata2_highsnr.best_values['exponent']
    specdata2_err_highsnr = resultpowdata2_highsnr.params['exponent'].stderr
    
    specfitdatageo_highsnr = -resultpowdatageo_highsnr.best_values['exponent']    
    specdatageo_err_highsnr = resultpowdatageo_highsnr.params['exponent'].stderr


    plt.figure(numfig5)
    plt.subplot(3,4,subpl_ind5)
    plt.plot(1000.*freqms_highsnr,taus2_order,'g*',markersize=9.0,linewidth=1.5,label=r'$\tau2$ from high snr data')
    plt.plot(1000.*freqms_highsnr, resultpowdata2_highsnr.best_fit, 'g-')
    plt.plot(1000.*freqms_highsnr,tausgeo,'r*',markersize=9.0,linewidth=1.5,label=r'$\tau_{geo}$')
    plt.plot(1000.*freqms_highsnr, resultpowdatageo_highsnr.best_fit, 'r-',label=r'$\alpha_{geo} = %.2f \pm %.2f$' %(specfitdatageo_highsnr,specdatageo_err_highsnr))  
    plt.title(r'$\alpha = %.2f \pm %.2f$' %(specfitdata2_highsnr,specfitdatageo_highsnr),fontsize=16)
    plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(r'$\nu$ MHz',fontsize=16)
    plt.ylabel(r'$\tau_2$ (sec)',fontsize=16)
    plt.legend(fontsize = 9)

##PLOT SIGMA##  #
#w50 = np.loadtxt(dri.skip_first_col('w50.list2'),skiprows=1)
w50s = np.array([ 12.2,   4.1,   6.2,   6.9,   5.4,  11. ,   4.8,   7.5,   4. ,14.8,   9. ,  14.8,  17.4])
freqw50s = np.hstack((0.408*np.ones(4),1.4*np.ones(3),0.408,1.4,0.408*np.ones(4)))

w50 = w50s[pulseind]/1000.
w50sig = w50/(2*np.sqrt(2*np.log(2)))
freqw50 = freqw50s[pulseind] 


if subpl_ind5 >= sp:
    subpl_ind6 = int(subpl_ind5-sp)+1
else:
    subpl_ind6 = int(subpl_ind5)+1

if npch + 6 == sp:    
    numfig6 = (npch+6)/sp
else:
    numfig6 = (npch+6)/sp + 1


figg = plt.figure(numfig6, figsize=(16,10))
figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05)    
plt.subplot(3,4,subpl_ind6)
plt.errorbar(1000.*freqms_highsnr,bestpT_highSNR[0]*pbs,yerr = bestpT_std_highSNR[0]*pbs, fmt = 'm*',markersize=9.0,capthick=2,linewidth=1.5,alpha=alfval)
#plt.plot(1000.*freqms_highsnr, sigmapow.best_fit, 'm-', label= r'$\beta$ = %.2f' %sigmapow.best_values['k'])
plt.plot(1000.*freqms_highsnr,powout.best_fit*pbs,'b-', alpha=alfval,label='pow = %.2f' %powout.best_values['exponent'])
#plt.plot(expout.best_fit*p2bs,1000.*freqms_highsnr, 'g-', label='exp = %.2f' %expout.best_values['decay'])
#plt.plot(1000.*freqms_highsnr, linout.best_fit, 'r-', label='m = %.2f' %linout.best_values['slope'])
plt.plot(1000.*freqms_highsnr,quadout.best_fit*pbs,'c-',alpha=alfval, label='a,b = %.2f,%.2f' %(quadout.best_values['a'],quadout.best_values['b']))
#plt.plot(220,w50sig,'k^',label=r'\$sigma$ at %d MHz' %(freqw50*1000))
plt.axhline(y=w50sig,color='k',label=r'$\sigma$ \rm{at} %d MHz' %(freqw50*1000))
plt.ylabel(r'$\sigma$ (sec)')
#plt.annotate(r'$\beta$ = %.2f' '\n m = %.2f' %(sigmapow.best_values['k'], linout.best_values['slope']),xy=(np.max(1000*freqms),np.max(bestpT[0])),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
#plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(bestpT[0])),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=11)
plt.xticks(fontsize=10)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
#plt.ylabel(r'$\sigma$',fontsize=16)
plt.legend(fontsize = 9, loc='best')


if subpl_ind6 >= sp:
    subpl_ind7 = int(subpl_ind6-sp)+1
else:
    subpl_ind7 = int(subpl_ind6)+1

if npch + 7 == sp:    
    numfig7 = (npch+7)/sp
else:
    numfig7 = (npch+7)/sp + 1

##PLOT DM##  

fig0 = plt.figure(numfig7, figsize=(16,10))
ax0 = fig0.add_subplot(3,4,subpl_ind7)

#fig0 = plt.figure(NF,figsize=(12,4))
#ax0 = plt.subplot(1,3,2)
#plt.errorbar(bestpT_highSNR[1]*pbs,1000.*freqms_highsnr, xerr = bestpT_std_highSNR[1]*pbs, fmt = 'b*',markersize=9.0,capthick=2,linewidth=1.5)
#plt.plot(1000.*freqms_highsnr, quadout_mu.best_fit*p2bs, 'b-', label='a,b = %.2f,%.2f' %(quadout_mu.best_values['a'],quadout_mu.best_values['b']))
#plt.plot(1000.*freqms_highsnr, powout_mu.best_fit*p2bs, 'm-', label='pow = %.2f' %powout_mu.best_values['exponent'])
#plt.plot(quadout_mu.best_fit*pbs,1000.*freqms_highsnr, 'b-', label='a,b = %.2f,%.2f' %(quadout_mu.best_values['a'],quadout_mu.best_values['b']))
#plt.plot(bestpT_highSNR[1]*pbs,1000*freqms_highsnr, 'k*-', label=r'DM: %.4f pc.cm^{-3}' %(DMval))
#plt.errorbar(delmuarray,delnuarray,xerr=delmu_stdarray, fmt='k*', label='DM rough check: %.4f' %DMcheck)
#plt.plot(DMmodelfit,delnuarray, 'm-', label=r'DM: %.4f pc.cm^{-3}' %(DMval))
#plt.errorbar(delmuarray,freqMHz_highsnr, xerr=delmu_stdarray, fmt='k*', label='DM rough check: %.4f' %DMcheck)
plt.errorbar(delmuarray,freqMHz_highsnr, xerr=delmu_stdarray, fmt=markr, alpha = alfval)
plt.plot(DMmodelfit,freqMHz_highsnr, 'k-', label=r'DM: $%.4f \pm %.4f$ $\rm{pc.cm}^{-3}$' %(DMval,DMvalstd), alpha = alfval)

#plt.plot(powout_mu.best_fit*p2bs,1000.*freqms_highsnr, 'm-', label='pow = %.2f' %powout_mu.best_values['exponent'])
#plt.xlabel(r'$\mu_H - \mu_{ch}$ (sec)', fontsize =18)
plt.xlabel(r'$\Delta \mu$ (sec)', fontsize =18)
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.title('%s' %pulsar)
#ax0.xaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
#plt.ylabel(r'$1/\nu_{ch}^2 - 1/\nu_H^2 $ MHz',fontsize=16)
plt.ylabel(r'$\nu$ (MHz)',fontsize=16)
plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
#plt.tight_layout()
plt.legend(fontsize = 12, loc='best')


delnuarray = [-(1/freqmsMHz[-1]**2-1/freqmsMHz[i]**2) for i in range(npch)]
delmuarray = [(bestpT_highSNR[1][-1] - bestpT_highSNR[1][i])*pbs for i in range(npch)]


## PLOT A ##

if subpl_ind7 >= sp:
    subpl_ind8 = int(subpl_ind7-sp)+1
else:
    subpl_ind8 = int(subpl_ind7)+1

if npch + 8 == sp:    
    numfig8 = (npch+8)/sp
elif npch + 8 == 2*sp:
    numfig8 = (npch+8)/sp
else:
    numfig8 = (npch+8)/sp + 1


plt.figure(numfig8,figsize=(16,10))
plt.subplot(3,4,subpl_ind8)
plt.errorbar(1000.*freqms_highsnr, bestpT_highSNR[2], yerr = bestpT_std_highSNR[2], fmt = 'g*',markersize=9.0,capthick=2,linewidth=1.5, alpha=alfval)
plt.title('Amplitude')
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
#plt.tight_layout()
plt.ylabel(r'$A$',fontsize=16)


## PLOT DC ##


if subpl_ind8 >= sp:
    subpl_ind9 = int(subpl_ind8-sp)+1
else:
    subpl_ind9 = int(subpl_ind8)+1

if npch + 9 == sp:    
    numfig9 = (npch+9)/sp
elif npch + 9 == 2*sp:
    numfig9 = (npch+9)/sp
else:
    numfig9 = (npch+9)/sp + 1


plt.figure(numfig9, figsize=(16,10))
plt.subplot(3,4,subpl_ind9)
plt.errorbar(1000.*freqms_highsnr, bestpT_highSNR[3], yerr = bestpT_std_highSNR[3], fmt = 'k*',markersize=9.0,capthick=2,linewidth=1.5,alpha=alfval)
plt.title('DC offset')
plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(bestpT_highSNR[3])),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlabel(r'$\nu$ MHz',fontsize=16)
#plt.tight_layout()
plt.ylabel(r'$DC$',fontsize=16)


## PLOT SNR AND TAU ERROR ###

if subpl_ind9 >= sp:
    subpl_ind10 = int(subpl_ind9-sp)+1
else:
    subpl_ind10 = int(subpl_ind9)+1

if npch + 10 == sp:    
    numfig10 = (npch+10)/sp
elif npch + 10 == 2*sp:
    numfig10 = (npch+10)/sp
else:
    numfig10 = (npch+10)/sp + 1


fig = plt.figure(numfig10, figsize=(16,10))
ax1 = fig.add_subplot(3,4,subpl_ind10)
ax2 = ax1.twinx()
ax1.axhline(y=SNRcutoff, xmin=0, xmax=len(freqmsMHz), linewidth=1.0, color = 'r',label='SNR cut = %.2f' %SNRcutoff,alpha=alfval)
ax2.axhline(y=tauerrorcutoff*100, xmin=0, xmax=len(freqmsMHz), linewidth=1.0, color = 'm',alpha=alfval,label=r'$\tau$' ' cut = %.1f' %(100*tauerrorcutoff))
ax1.plot(freqmsMHz, comp_SNRs,'g^',markersize=7, label='SNR',alpha=alfval)
ax2.plot(freqmsMHz,np.array(lmfittausstds)/np.array(obtainedtaus)*100,'b*',markersize=9, label=r'$\tau$',alpha=alfval)
plt.title(r'SNR and $\tau$ error ($\%$)')
ax1.set_xlabel(r'$\nu$ MHz',fontsize=16)
#ax1.set_ylabel('SNR',fontsize=16)
#ax2.set_ylabel(r'$\tau$' ' error (%)',fontsize=16)
ax2.set_yscale('log')
leg1 = ax1.legend(loc='lower left',fontsize=9)
leg2 = ax2.legend(loc='upper right',fontsize=9)
leg1.get_frame().set_alpha(0.7)  #make legends transparent
leg2.get_frame().set_alpha(0.7)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
#plt.tight_layout()


##plt.figure(10)
##plt.plot(0.1*np.sum(profiles,axis=0), 'k', lw=3.0)
##plt.axhline(y=0.1*np.max(np.sum(profiles,axis=0)), label='Max = %.2f' %np.max(np.sum(profiles,axis=0)))
##[plt.plot(profiles[i]) for i in range(npch)]
##plt.title('Shifted unscattered (Sum/10)')
##plt.xlim(xmax=2.5*np.argmax(np.sum(profiles,axis=0)))
##plt.legend()
##plotnametau = '%s_Sum_shifted_unscattered.png'  % (pulsar)
##picpathtau = newpath
##fileoutputtau = os.path.join(picpathtau,plotnametau)
##plt.savefig(fileoutputtau, dpi=150)
##plt.close(10)
##
##plt.figure(11)
##plt.plot(np.sum(datas, axis=0), 'r', lw=3.0)
##plt.axhline(y=np.max(np.sum(datas,axis=0)), label='Max = %.2f' %np.max(np.sum(datas,axis=0)))
##plt.title('Shifted summed data')
##plt.legend()
##plotnametau = '%s_Sum_shifted_data.png'  % (pulsar)
##picpathtau = newpath
##fileoutputtau = os.path.join(picpathtau,plotnametau)
##plt.savefig(fileoutputtau, dpi=150)
##plt.close(11)
#
#### PLOT ARTICLE TYPE FIGURE, 1 by 3
#
#NF = 10
#
#figof3 = plt.figure(NF,figsize=(12,4))
#figof3.subplots_adjust(left = 0.065, right = 0.97, wspace=0.25, bottom=0.15)
#plt.subplot(1,3,3)
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.plot(freqMHz_highsnr,fluxes_highsnr,markr,alpha=alfval,markersize=11)
#plt.plot(freqMHz_highsnr, fluxes_highsnr+climb_highsnr,prof,linewidth=1.0)
##plt.plot(freqMHz_highsnr, unscatflux,'go',markersize=10.0, alpha=0.5)
#plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=alfval2, facecolor=lcol)
#plt.title('%s' %(pulsar))
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)
#plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
#plt.ylabel(r'Calibrated flux (mJy)',fontsize=16)
#
#plt.subplot(1,3,1)
#plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,alpha=alfval, markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
#plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-')
#plt.title('%s' %pulsar)
#plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
#plt.ylabel(r'$\tau$ (sec)',fontsize=16)
#plt.legend(fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(xmin=(freqms[0])*950,xmax=(freqms[-1])*1050)
#plt.yticks(fontsize=12)
#plt.xticks(ticksMHz,ticksMHz,fontsize=12)
##plt.tight_layout()
#
#fig0 = plt.figure(NF,figsize=(12,4))
#ax0 = plt.subplot(1,3,2)
#plt.errorbar(delmuarray,freqMHz_highsnr, xerr=delmu_stdarray, fmt=markr,alpha=alfval)
#plt.plot(DMmodelfit,freqMHz_highsnr, 'm-', label=r'$\Delta$DM: $%.4f \pm %.4f$ $\rm{pc.cm}^{-3}$' %(DMval,DMvalstd))
#plt.xlabel(r'$\Delta \mu$ (sec)', fontsize =18)
#plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
#plt.yticks(fontsize=12)
#plt.xticks(fontsize=12)
#plt.title('PSR %s' %pulsar)
#plt.ylabel(r'$\nu$ (MHz)',fontsize=16)
#plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
##plt.tight_layout()
#plt.legend(fontsize = 12, loc='best')


NF = 11

nrr = 2

figof3 = plt.figure(NF,figsize=(12,4))
figof3.subplots_adjust(left = 0.075, right = 0.97, wspace=0.25, bottom=0.15)
plt.subplot(1,3,1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(profilexaxis,data_highsnr[nrr],'k',alpha = 0.1 )
plt.plot(profilexaxis,model_highsnr[nrr],prof, lw=1.5, label =r'$\tau: %.0f \pm %.0f$ ms' %(1000*taus_highsnr[nrr]*pulseperiod/nbins, 1000*lmfitstds_highsnr[nrr]*pulseperiod/nbins))
plt.title('PSR %s at %.1f MHz' %(pulsars[nrr], freqMHz_highsnr[nrr]))
#plt.annotate(r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[4]*pulseperiod/nbins, lmfitstds_highsnr[4]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[4])),xycoords='data',xytext=(0.4,0.9),textcoords='axes fraction',fontsize=12)
plt.ylim(ymax=1.35*np.max(data_highsnr[nrr]))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time (sec)',fontsize=16)
plt.ylabel('intensity (mJy)',fontsize=16)
plt.xlim(xmin=0,xmax=pulseperiod)
plt.legend(fontsize=14)


#freq2KS = np.arange(1000.*freqms_highsnr[0],327.0,10)

def round_down(num, divisor):
    return num - (num%divisor)

ax = plt.subplot(1,3,2)
stf,endf= 0,0

if meth in 'iso':
    if pulsar == 'J0614+2229':
        stf,endf = 100, 350
    if pulsar == 'J0742-2822':
        stf,endf = round_down(freqMHz_highsnr[0],10), round(freqMHz_highsnr[-1],10)
    if pulsar == 'J1909+1102':
        stf,endf = round_down(freqMHz_highsnr[0],10), round(freqMHz_highsnr[-1],10)
    if pulsar == 'J1917+1353':
        stf,endf = 100, 200
    if pulsar == 'J1922+2110':
        stf,endf = 100, 250
    if pulsar == 'J1935+1616':
        stf,endf = 100, 250
    if pulsar == 'J2305+3100':
        stf,endf = 100, 200
    if pulsar in ('J0040+5716, J0543+2329, J1851+1259, J1913-0440, J2257+5909'):
        stf,endf = freqMHz_highsnr[0], freqMHz_highsnr[-1]

if pulsar == 'J0614+2229' and meth in 'onedim':
   stf,endf = 100, 350
   plt.errorbar(327.0, 1.74/1000., yerr = 0.03/1000.,fmt='b^', label='K15', markersize=11.0)
   plt.errorbar(327.0, 0.0028, yerr=0.02/1000.,fmt='r*', label='Our fit to K15', markersize=12.0)
   plt.errorbar(111.,40./1000.,yerr=10/1000.,fmt='co',label='Kuzmin 2007',markersize=9.0)
   plt.ylim(ymax=1.7)
if pulsar == 'J0742-2822' and meth in 'onedim':
    stf,endf = freqMHz_highsnr[0], freqMHz_highsnr[-1]
    plt.errorbar(160, 24.5/1000.,yerr=2.8/1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
if pulsar == 'J1909+1102' and meth in 'onedim':
    stf,endf = freqMHz_highsnr[0], freqMHz_highsnr[-1]
    plt.errorbar(160, 27.0/1000.,yerr=7.0/1000., fmt = 'g^', label='Slee 1980', markersize=12.0)
    plt.errorbar(160, 26.5/1000.,yerr=8.1/1000., fmt = 'ro', label='Alukar 1986', markersize=7.0)
if pulsar == 'J1917+1353' and meth in 'onedim':
    stf,endf = 100, 200
    plt.errorbar(160, 11.7/1000.,yerr=1.9/1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
    plt.errorbar(160, 12./1000.,yerr=3./1000., fmt = 'ro', label='Slee 1980', markersize=7.0)
    plt.errorbar(102, 40./1000.,yerr=20./1000.,fmt= 'co', label='Kuzmin 2007', markersize=9.0)
if pulsar == 'J1922+2110' and meth in 'onedim':
    stf,endf = 100, 250
    plt.plot(102.75, 26/1000., 'm^', label='L15', markersize=12.0)    
    plt.errorbar(160, 96.8/1000.,yerr=50./1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
    plt.errorbar(243, 4.4/1000.,yerr=1.5/1000., fmt = 'k^', label=r'L\"{o}hmer 2004', markersize=12.0)
    plt.ylim(10**-3,0.5)
if pulsar == 'J1935+1616' and meth in 'onedim':
    stf,endf = 100, 250
    plt.errorbar(160, 21.7/1000.,yerr=1.6/1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
    plt.errorbar(160, 25./1000.,yerr=4./1000., fmt = 'ro', label='Slee 1980', markersize=7.0)
    plt.errorbar(102, 50./1000.,yerr=15./1000.,fmt= 'co', label='Kuzmin 2007', markersize=9.0)
    plt.errorbar(243, 4.6/1000.,yerr=0.2/1000., fmt = 'k^', label=r'L\"{o}hmer 2004', markersize=12.0)
    plt.ylim(2.5*10**-3,1.5)
if pulsar == 'J2305+3100' and meth in 'onedim':
    stf,endf = 100, 200
    plt.errorbar(160, 9.9/1000.,yerr=3.6/1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
    plt.errorbar(102, 13./1000.,yerr=3./1000.,fmt= 'co', label='Kuzmin 2007', markersize=9.0)
if pulsar in ('J0040+5716, J0543+2329, J1851+1259, J1913-0440, J2257+5909') and meth in 'onedim':
    stf,endf = freqMHz_highsnr[0], freqMHz_highsnr[-1]


plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,alpha=alfval, markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
#plt.plot(1000.*freqms_highsnr, powouttau_highsnr.best_fit, 'k-')
freq2KS = np.arange(float(stf),float(endf),0.10)

ticksMHzL = range(int(stf),int(endf)+20,20)
Amp = powouttau_highsnr.best_values['amplitude']
Expo = powouttau_highsnr.best_values['exponent']
plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), 'k-',alpha=alfval)
plt.title('PSR %s' %pulsar)
plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
plt.ylabel(r'$\tau$ (sec)',fontsize=16)
#plt.legend(fontsize=14,loc='best',numpoints=1)
if meth in 'onedim':
    handles, labels = ax.get_legend_handles_labels()
    labels.insert(1,labels.pop(-1))
    handles.insert(1,handles.pop(-1))
    if len(labels) > 3:
        fs = 13
    else:
        fs=14
    ax.legend(handles, labels, loc='best',numpoints=1, fontsize=fs)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.9*stf,1.1*endf)
#plt.xlim(xmin=(freqms[0])*950,xmax=(freqms[-1])*1050)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xticks(ticksMHzL,ticksMHzL,fontsize=12)
#plt.tight_layout()

if pulsar in ('J1917+1353','J1922+2110', 'J1935+1616'):
    fig0 = plt.figure(NF,figsize=(12,4))
    ax0 = plt.subplot(1,3,3)
    plt.errorbar(delmuarray,freqMHz_highsnr, xerr=delmu_stdarray, fmt=markr,alpha=alfval, markersize=9.0,label=r'$\Delta$DM:$%.4f \pm %.4f$ $\rm{pc.cm}^{-3}$' %(DMval,DMvalstd))
    plt.plot(DMmodelfit,freqMHz_highsnr, 'k-', alpha=alfval)
    plt.xlabel(r'$\Delta \mu$ (sec)', fontsize =18)
    plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.title('PSR %s' %pulsar)
    plt.ylabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
    #plt.tight_layout()
    plt.legend(fontsize = 12, loc='best',numpoints=1)
#
#
#def add_arrow(x1,col,ls):
#    """Add a vertical force arrow at x1 (in data coordinates)."""
#    plt.annotate('', xy=(freqMHz_highsnr[x1], fluxes_highsnr[x1]), xytext=(freqMHz_highsnr[x1], fluxes_highsnr[x1]+climb_highsnr[x1]), textcoords='data', arrowprops=dict(arrowstyle='<|-', color=col, linestyle=ls))
#
#
    
else:
    plt.subplot(1,3,3)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #plt.plot(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,'go',markersize=12)
    #plt.plot(freqMHz_highsnr, unscatflux,'co',markersize=7.0)
    plt.plot(freqMHz_highsnr, fluxes_highsnr+climb_highsnr,prof,linewidth=2.0)
    plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=alfval2,facecolor=lcol)
    eb = plt.errorbar(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,yerr=sigmaFlux, fmt=lcol,markersize=0.0, alpha=1.0,capthick=2,linewidth=1.5)
    eb[-1][0].set_linestyle(ls)
    plt.plot(freqMHz_highsnr,fluxes_highsnr,markr,alpha=alfval,markersize=11)
    plt.title('PSR %s' %(pulsar))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ylabel(r'mean flux density (mJy)',fontsize=16)
    plt.ylim(ymin=0.9*np.min(fluxes_highsnr+climb_highsnr-sigmaFlux),ymax=1.1*np.max(fluxes_highsnr+climb_highsnr+sigmaFlux))
    plt.xlim(0.98*freqMHz_highsnr[0],1.02*freqMHz_highsnr[-1])
##if meth in 'iso':
##for i in range(npch):
##    #    plt.axvline(freqMHz_highsnr[i])
##        add_arrow(i,lcol,ls)
#
#
#



"""Appendix Plots"""


figg = plt.figure(numFig,figsize=(16,10))
figg.subplots_adjust(left = 0.055, right = 0.98, wspace=0.35,hspace=0.35,bottom=0.05) 
    
if nch/8 == 1:
    bb = 1
elif nch == 4:
    bb = 1
else:
    bb = 2
print "bb = %d" %bb    

NF=25

#if pulsar == 'J1851+1259'  and meth == 'iso':
#    data_highsnr = np.delete(data_highsnr, (8), axis=0)
#    model_highsnr = np.delete(model_highsnr, (8), axis=0)
#    freqMHz_highsnr = np.delete(freqMHz_highsnr, (8), axis=0)
#    lmfitstds_highsnr = np.delete(lmfitstds_highsnr, (8), axis=0)
#    taus_highsnr = np.delete(taus_highsnr, (8), axis=0)

figof8 = plt.figure(NF,figsize=(16,8))
figof8.subplots_adjust(left = 0.06, right = 0.98, wspace=0.25, bottom=0.10, hspace=0.35)
count = 1
for i in range(0,npch,bb):
#for i in (0,2,4,6,9,11):   ##use for J1851 iso
    plt.subplot(2,4,count)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(profilexaxis,data_highsnr[i],'k',alpha = 0.1 )
    plt.plot(profilexaxis,model_highsnr[i],prof, lw=1.5, label =r'$\tau: %.0f \pm %.0f$ ms' %(taus_highsnr[i]*pulseperiod/nbins*1000, lmfitstds_highsnr[i]*pulseperiod/nbins*1000))
    plt.title('PSR %s at %.1f MHz' %(pulsars[i], freqMHz_highsnr[i]))
    #plt.annotate(r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[4]*pulseperiod/nbins, lmfitstds_highsnr[4]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[4])),xycoords='data',xytext=(0.4,0.9),textcoords='axes fraction',fontsize=12)
#    plt.ylim(ymin=-1000 + (-500*count),ymax=1.35*np.max(data_highsnr[i]))
    lowl = np.min(data_highsnr[0])/np.max(data_highsnr[0])
    plt.ylim(ymin=lowl*np.max(data_highsnr[i]),ymax=1.6*np.max(data_highsnr[i]))
#    plt.xticks([0,0.04,0.08,0.12,0.16],[0,0.04,0.08,0.12,0.16],fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
#    plt.xlim(0.2,0.7)
#    plt.xlim(0.6,1.0)
#    plt.xlim(0.3,0.9)
    plt.xlim(xmax=pulseperiod)
    plt.legend(fontsize=14)
    count += 1

plt.subplot(2,4,1)
plt.ylabel('intensity (mJy)',fontsize=16)
plt.subplot(2,4,5)
plt.ylabel('intensity (mJy)',fontsize=16)

for r in range(5,9):
    plt.subplot(2,4,r)
    plt.xlabel('time (sec)',fontsize=16)
    

#ax1 = plt.subplot(2,4,8)
#ax1.axis('off')


    

#tk = np.arange(0, 0.20, 0.04)
#
#count=1
#NF=25
#figof8 = plt.figure(NF,figsize=(16,4))
#figof8.subplots_adjust(left = 0.06, right = 0.98, wspace=0.25, bottom=0.12, hspace=0.35)
#count = 1
#for i in range(0,npch,1):
#    plt.subplot(1,4,count)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    plt.plot(profilexaxis,data_highsnr[i],'y',alpha = 0.25 )
#    plt.plot(profilexaxis,model_highsnr[i],prof, alpha = 0.5, label =r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[i]*pulseperiod/nbins, lmfitstds_highsnr[i]*pulseperiod/nbins))
#    plt.title('%s at %.1f MHz' %(pulsars[i], freqMHz_highsnr[i]))
#    #plt.annotate(r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[4]*pulseperiod/nbins, lmfitstds_highsnr[4]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[4])),xycoords='data',xytext=(0.4,0.9),textcoords='axes fraction',fontsize=12)
##    plt.ylim(ymin=-1000 + (-500*count),ymax=1.35*np.max(data_highsnr[i]))
#    lowl = np.min(data_highsnr[0])/np.max(data_highsnr[0])
#    plt.ylim(ymin=lowl*np.max(data_highsnr[i]),ymax=1.6*np.max(data_highsnr[i]))
#    plt.xticks(tk,tk,fontsize=14)
#    plt.yticks(fontsize=14)
##    plt.xlim(0.2,0.7)
#    plt.xlim(xmax=pulseperiod)
#    plt.legend(fontsize=11)
#    plt.xlabel('time (s)',fontsize=16)
#    count += 1
#
#plt.subplot(1,4,1)
#plt.ylabel('intensity (mJy)',fontsize=16)


###PLOT FOR PSR J1913-0440

if pulsar == 'J1913-0440':

    NF = 13
    
    figof3 = plt.figure(NF,figsize=(12,9))
    figof3.subplots_adjust(left = 0.075, right = 0.97, wspace=0.35, bottom=0.12, hspace=0.35)
    count = 1
    for i in (0,5,12):
        plt.subplot(2,3,count)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if datac == 'comm':
            plt.plot(profilexaxis,data_highsnr[i],'k',alpha = 0)
            if meth == 'iso':
                plt.plot(profilexaxis,model_highsnr[i],'k:', alpha = 1.0, lw=1.0)
        else:
            plt.plot(profilexaxis,data_highsnr[i],'k',alpha = 0.1 )
            plt.plot(profilexaxis,model_highsnr[i],prof, alpha = 1.0, lw=1.5,label =r'$\tau: %.0f \pm %.0f$ ms' %(1000*taus_highsnr[i]*pulseperiod/nbins, 1000*lmfitstds_highsnr[i]*pulseperiod/nbins))
        plt.title('PSR %s at %.1f MHz' %(pulsars[i], freqMHz_highsnr[i]))
        #plt.annotate(r'$\tau: %.4f \pm %.4f$ sec' %(taus_highsnr[4]*pulseperiod/nbins, lmfitstds_highsnr[4]*pulseperiod/nbins),xy=(np.max(profilexaxis),np.max(data_highsnr[4])),xycoords='data',xytext=(0.4,0.9),textcoords='axes fraction',fontsize=12)
    #    plt.ylim(ymin=-1000 + (-500*count),ymax=1.35*np.max(data_highsnr[i]))
        ymaxx = [12000,0,0,0,0,30000,0,0,0,0,0,0,50000]        
        plt.ylim(ymin=-0.2*np.max(data_highsnr[i]),ymax=ymaxx[i])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('time (sec)',fontsize=16)
        plt.xlim(0.3,0.6)
        plt.legend(fontsize=14)
        count += 1
    tix = np.arange(110,210,20)
    plt.subplot(2,3,1)
    plt.ylabel('intensity (mJy)',fontsize=16)
    
    plt.subplot(2,3,4)
    freq2KS = np.arange(float(100),float(200),0.10)
    Amp = powouttau_highsnr.best_values['amplitude']
    Expo = powouttau_highsnr.best_values['exponent']
    if datac == 'comm':
            if meth == 'iso':
                mm = 'k:'
            else:
                mm = 'k-.'
#            plt.plot(1000.*freqms_highsnr,taussec_highsnr,mm,alpha=alfval, markersize=9.0,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
            plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), mm)
    else:
        plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,alpha=alfval, markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
        plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), 'k-')
        if meth in 'onedim':    
                plt.errorbar(160, 0.0167,yerr=1.8/1000., fmt = 'g^', label='Alukar 1986', markersize=12.0)
                plt.errorbar(160, 32./1000.,yerr=5./1000., fmt = 'ro', label='Slee 1980', markersize=7.0)
                plt.errorbar(102, 35./1000.,yerr=15./1000.,fmt= 'co', label='Kuzmin 2007', markersize=9.0)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ylabel(r'$\tau$ (sec)',fontsize=16)
    plt.legend(fontsize=12,numpoints=1,loc='best')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(xmin=95,xmax=(freqms[-1])*1050)
    plt.yticks(fontsize=12)
    plt.xticks(tix,tix,fontsize=12)
    plt.ylim(ymin=5*10**-4)
    
    
    
    #plt.tight_layout()
    
    fig0 = plt.figure(NF,figsize=(12,4))
    ax0 = plt.subplot(2,3,5)
    if datac == 'comm': 
        if meth == 'iso':
                mm = 'k:'
        else:
                mm = 'k-.'
        plt.plot(DMmodelfit,freqMHz_highsnr,mm, alpha=alfval) 
    else:
        plt.errorbar(delmuarray,freqMHz_highsnr, xerr=delmu_stdarray, fmt=markr,alpha=alfval, markersize=9.0,label=r"$\Delta$DM:$%.4f \pm %.4f$ $\rm{pc.cm}^{-3}$" %(DMval,DMvalstd))
        plt.plot(DMmodelfit,freqMHz_highsnr, 'k-', alpha=alfval)    
    plt.xlabel(r'$\Delta \mu$ (sec)', fontsize =18)
    plt.annotate('%s' %pulsar,xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    #plt.title('%s' %pulsar)
    plt.ylabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
    #plt.tight_layout()
    plt.legend(fontsize = 12, loc='best',numpoints=1)
    
    
    
    
    plt.subplot(2,3,6)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif') 
    if datac == 'comm':
        if meth == 'iso':
                mm = 'k:'
        else:
                mm = 'k-.'
        plt.plot(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,mm,)    
    else:
        plt.plot(freqMHz_highsnr, fluxes_highsnr+climb_highsnr,prof,linewidth=2.0)
        plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=alfval2,facecolor=lcol)
        eb = plt.errorbar(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,yerr=sigmaFlux, fmt=lcol,markersize=0.0, alpha=1.0,capthick=2,linewidth=1.5)
        eb[-1][0].set_linestyle(ls)
        plt.plot(freqMHz_highsnr,fluxes_highsnr,markr,alpha=alfval,markersize=11)
        #plt.plot(freqMHz_highsnr, fluxes_highsnr+bestpT_highSNR[3],'r-',linewidth=1.0,label='DC added')
    #plt.title('%s' %(pulsar))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(110,190)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ylabel(r'Calibrated flux (mJy)',fontsize=16)
#    for i in range(npch):
#    #    plt.axvline(freqMHz_highsnr[i])
#        add_arrow(i,lcol,ls)
##    plt.legend(fontsize = 12, loc='best')
    
    
    

###PLOT FOR J0117+5914

if pulsar in ('J0614+2229','J0117+5914'):
    
    #Subplot1

    NF = 14
    
    farr = 1 

    figof3 = plt.figure(NF,figsize=(12,4))
    figof3.subplots_adjust(left = 0.075, right = 0.97, wspace=0.25, bottom=0.15)
    plt.subplot(1,3,1)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    if datac == 'comm':
        plt.plot(profilexaxis,data_highsnr[farr],'k',alpha = 0)
        if meth == 'iso':
            plt.plot(profilexaxis,model_highsnr[farr],'k:', alpha = 1.0, lw=1.0)
    else:
        plt.plot(profilexaxis,data_highsnr[farr],'k',alpha = 0.1 )
        plt.plot(profilexaxis,model_highsnr[farr],prof, alpha = 1.0, lw=1.5,label =r'$\tau: %.0f \pm %.0f$ ms' %(1000*taus_highsnr[farr]*pulseperiod/nbins, 1000*lmfitstds_highsnr[farr]*pulseperiod/nbins))
    #plt.plot(profilexaxis,data_highsnr[2],'k',alpha = 0.1 )
    #plt.plot(profilexaxis,model_highsnr[2],prof, lw=1.5, label =r'$\tau: %.1f \pm %.1f$ ms' %(1000*taus_highsnr[2]*pulseperiod/nbins, 1000*lmfitstds_highsnr[2]*pulseperiod/nbins))
    plt.title('PSR %s at %.1f MHz' %(pulsars[farr], freqMHz_highsnr[farr]))
    plt.ylim(-100,800)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('time (sec)',fontsize=16)
    plt.ylabel('intensity (mJy)',fontsize=16)
    plt.xlim(xmin=0,xmax=pulseperiod)
    plt.legend(fontsize=14)

    #Subplot2

    ax = plt.subplot(1,3,2)
    if pulsar == 'J0614+2229':
        stf,endf = 100, 350
    else:
        stf,endf = 100, 200
    freq2KS = np.arange(float(stf),float(endf),0.10)
    Amp = powouttau_highsnr.best_values['amplitude']
    Expo = powouttau_highsnr.best_values['exponent']
    if datac == 'comm':
            if meth == 'iso':
                mm = 'k:'
            else:
                mm = 'k-.'
#            plt.plot(1000.*freqms_highsnr,taussec_highsnr,mm,alpha=alfval, markersize=9.0,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
            plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), mm,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
    else:
#        plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), 'k-',alpha=alfval)
        plt.errorbar(1000.*freqms_highsnr,taussec_highsnr,yerr=lmfitstds_highsnr*pulseperiod/nbins,fmt=markr,alpha=alfval, markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.1f \pm %.1f$' %(specfitdata_highsnr,specdata_err_highsnr))
        
        plt.plot(freq2KS, Amp*np.power(freq2KS/1000.,Expo), 'k-')
    ticksMHzL = range(110,endf+10,20)
    plt.title('PSR %s' %pulsar)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ylabel(r'$\tau$ (sec)',fontsize=16)
    plt.legend(fontsize=14,loc='best',numpoints=1)
    if pulsar == 'J0614+2229' and meth in 'onedim' and datac in 'comm':
        plt.errorbar(327.0, 1.74/1000., yerr = 0.03/1000.,fmt='b^', label='K15', markersize=11.0)
        plt.errorbar(327.0, 0.0028, yerr=0.02/1000.,fmt='r*', label='Our fit to K15', markersize=12.0)
        plt.errorbar(111.,40./1000.,yerr=10/1000.,fmt='co',label='Kuzmin 2007',markersize=9.0)
        plt.ylim(ymax=1.7)  
    if meth in 'onedim' and datac in 'comm':
        ax = plt.subplot(1,3,2)
        handles, labels = ax.get_legend_handles_labels()
        labels.insert(0,labels.pop(-2))
        handles.insert(0,handles.pop(-2))        
        labels.insert(1,labels.pop(-1))
        handles.insert(1,handles.pop(-1))
        ax.legend(handles, labels, loc='best',numpoints=1, fontsize=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.9*stf,1.1*endf)
    plt.yticks(fontsize=12)
    plt.xticks(fontsize=12)
    plt.xticks(ticksMHzL,ticksMHzL,fontsize=12)
    #plt.tight_layout()




    #Subplot3

    plt.subplot(1,3,3)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    if datac == 'comm':
        if meth == 'iso':
                mm = 'k:'
        else:
                mm = 'k-.'
        plt.plot(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,mm,lw=2.0)
        plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=0.3,facecolor='k')
        eb = plt.errorbar(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,yerr=sigmaFlux, fmt='k*',markersize=0.0, alpha=1.0,capthick=2,linewidth=1.5)
        eb[-1][0].set_linestyle(ls)
    else:
        plt.plot(freqMHz_highsnr, fluxes_highsnr+climb_highsnr,prof,linewidth=2.0)
        plt.fill_between(freqMHz_highsnr,fluxes_highsnr,fluxes_highsnr+climb_highsnr,alpha=alfval2,facecolor=lcol)
        eb = plt.errorbar(freqMHz_highsnr,fluxes_highsnr+climb_highsnr,yerr=sigmaFlux, fmt=lcol,markersize=0.0, alpha=1.0,capthick=2,linewidth=1.5)
        eb[-1][0].set_linestyle(ls)
        plt.plot(freqMHz_highsnr,fluxes_highsnr,markr,alpha=alfval,markersize=11)    
    plt.title('PSR %s' %(pulsar))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
    plt.ylabel(r'mean flux density (mJy)',fontsize=16)
    plt.ylim(ymin=0.9*np.min(fluxes_highsnr+climb_highsnr-sigmaFlux),ymax=150)
    plt.xlim(110,190)










### CONFIDENCE INTERVALS

#levels = [0.997, 0.95, 0.68, 0.00]
#
#for k in range(1):
#    plt.figure(500*(k+1),figsize=(8,8))
##    limitx = (taussec_highsnr[k]/pbs - lmfitstds_highsnr[k],taussec_highsnr[k]/pbs + 8*lmfitstds_highsnr[k])             
##    limity = (bestpT_highSNR[2][k] - bestpT_std_highSNR[2][k], bestpT_highSNR[2][k] + 8*bestpT_std_highSNR[2][k])
#    limitx = (0.05/pbs,0.20/pbs)
##    limity = (1200,3200)    
#    limity = (0.002/pbs,0.003/pbs)
#    if meth in 'onedim':     
#        cx, cy, grid = lmfit.conf_interval2d(results[k], 'tau1','sigma',100,100,limits=(limitx,limity))
#    else:
#        cx, cy, grid = lmfit.conf_interval2d(results[k], 'tau','A',100,100,limits=(limitx,limity))
##    CS = plt.contour(cx*pbs, cy, grid, levels)
#    CS = plt.contour(cx*pbs, cy*pbs, grid, levels)
#    plt.clabel(CS, fmt='%.2f', colors='k', fontsize=14)    
#    plt.xlabel('tau (sec)')
##    plt.ylabel('A')
#    plt.ylabel('sigma')
#

##### Plot the residuals and autocorrelations of those residuals
##### PLOT RESIDUALS

#for i in range(npch/sp +1):
#    plt.figure(numfig10+(i+1), figsize=(16,10))
#    
#for j in range(npch):
#    if j+1 == sp:
#        numFig = (j+1)/sp
#        subplotcount = sp
#    else: 
#        numFig = (j+1)/sp + 1
#        if j+1 < sp:
#            subplotcount = j+1
#        else: 
#            subplotcount = j+1 - sp
#    figres1 = plt.figure(numfig10 + numFig)
#    figres1.subplots_adjust(left = 0.055, right = 0.98, wspace=0.25,hspace=0.35,bottom=0.05)          
#    plt.subplot(3,4,subplotcount)
##    plt.tight_layout()
##    plt.plot(profilexaxis, data_highsnr[j] - model_highsnr[j],'b', alpha=0.5, label = 'residuals')
#    plt.hist(data_highsnr[j] - model_highsnr[j],alpha=alfval)
##    plt.plot(profilexaxis, psr.autocorr(data_highsnr[j] - model_highsnr[j]),'r', label = 'ACF')
#    plt.title('%s at %.1f MHz' %(pulsars[j], freqMHz_highsnr[j]))
##    plt.annotate(r'residuals',xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
#    plt.xticks(fontsize=12)
#    plt.yticks(fontsize=12)
##    plt.xlim(-0.1*pulseperiod,plt.xlim()[1])
#    plt.xlabel('intensity residuals')

### PLOT ACFs

#if npch!=sp:
#    for i in range(npch/sp+1):
#        plt.figure(numfig10+(npch/sp+1)+(i+1), figsize=(16,10))
#    
#for j in range(npch):
#    if j+1 == sp:
#        numFig = (j+1)/sp
#        subplotcount = sp
#        if npch == sp:
#            numFig = numFig-1
#    else: 
#        numFig = (j+1)/sp + 1
#        if npch == sp:
#            numFig = numFig-1
#        if j+1 < sp:
#            subplotcount = j+1
#        else: 
#            subplotcount = j+1 - sp
#    figres = plt.figure(numfig10 +(npch/sp+1) +numFig, figsize=(16,10))         
#    figres.subplots_adjust(left = 0.055, right = 0.98, wspace=0.25,hspace=0.35,bottom=0.05)        
#    plt.subplot(3,4,subplotcount)    
##    plt.tight_layout()
##    plt.plot(profilexaxis, data_highsnr[j] - model_highsnr[j],'b', alpha=0.5, label = 'residuals')
#    plt.plot(profilexaxis[1:], psr.autocorr(data_highsnr[j] - model_highsnr[j])[1:],'r', label = 'ACF',alpha=alfval)
#    plt.title('%s at %.1f MHz' %(pulsars[j], freqMHz_highsnr[j]))
##    plt.annotate(r'residuals',xy=(np.max(1000*freqms),np.max(obtainedtausec)),xycoords='data',xytext=(0.5,0.7),textcoords='axes fraction',fontsize=14)
#    plt.xticks(fontsize=12)
#    plt.yticks(fontsize=12)
#    plt.legend()
#    plt.xlabel('sec')
#    plt.xlim(-0.1*pulseperiod,plt.xlim()[1])
#
## PLOT TAU REFERENCED TO KRISHNAKUMAR

#Krishna_values_ms = np.array([1.74, 0.71, 1.35, 0.19, 0.36, 2.3, 3.21])
#Krishna_errors_ms = np.array([0.03, 0.01, 0.02, 0.01, 0.01, 0.1, 0.02])
#
#Krishna_pulsars = ['J0614+2229', 'J0742-2822','J1909+1102','J1913-0440','J1917+1353','J1922+2110','J1935+1616']
#
#if pulsar in Krishna_pulsars:
#	Kh_ind = Krishna_pulsars.index(pulsar)
#
#	Krishna = Krishna_values_ms[Kh_ind]
#	Krishna_Err = Krishna_errors_ms[Kh_ind]
#
#	taus327 = []
#	for i in range(npch):
#	    tau327 = psr.tauatfreq(1000.*freqms_highsnr[i],taussec_highsnr[i],327.0,specfitdata_highsnr)
#	    taus327.append(tau327)
#
#	taus327 = np.array(taus327)
#
#	#mean_err = (np.max(taus327)-np.min(taus327))/(2*np.sqrt(npch))  #Def for error in average value used: ave_std = range/(2*sqrt(N)), with range = max_value - min_value, and N number of data points.
#
#	spec_link = psr.linking_specind(taussec_highsnr[-1],1000.*freqms_highsnr[-1],Krishna/1000.,327)
#	spec_link_mean = psr.linking_specind(np.mean(taussec_highsnr),1000.*np.mean(freqms_highsnr),Krishna/1000.,327)
#
#
#	plt.figure(12)
#	plt.plot(1000.*freqms_highsnr,1000.*taus327,markr,alpha=alfval,markersize=11.0,label=r'LOFAR extrapolation from each $\nu_{obs}$')
#	plt.axhline(y=Krishna,lw=2.0,label='Krishnakumar et al.: %.2f ms' %(Krishna))
#	plt.fill_between(1000.*freqms_highsnr,Krishna-Krishna_Err, Krishna+Krishna_Err,facecolor='b',alpha=0.4 )
#	plt.axhline(y=1000*np.mean(taus327),color='r',label='LOFAR mean value at 327 MHz: %.2f ms' % (1000*np.mean(taus327)))
#	#plt.fill_between(1000.*freqms_highsnr,1000*np.mean(taus327)-1000*mean_err,1000*np.mean(taus327)+1000*mean_err)
#	plt.ylabel(r'$\tau$ at 327 MHz (ms)')
#	plt.xlabel(r'$\nu_{obs}$ MHz',fontsize=16)
#	plt.text(1.1*plt.xlim()[0],textpos*plt.ylim()[1],r'%s, link high $\nu$ %s: $\alpha$ = %.2f, link ave. $\nu$: $\alpha$ = %.2f' %(pulsar,meth,spec_link,spec_link_mean)) 
#	plt.legend(loc='best', fontsize=9)
#
#if temp is not None:
#  bf, A, S = psr.FitTemplate(tempdata,nbins)
#
#  plt.figure()
#  plt.plot(tempdata,'y')
#  plt.plot(bf)
#
##spec = np.zeros(10)
##for i in range(10):
##    bb = 11 - (i+2)
##    ee = 11
##    print freqms[bb:ee]
##    powparstau = powmod.guess(obtainedtausec[bb:ee],x=freqms[bb:ee])
##    powouttau = powmod.fit(obtainedtausec[bb:ee],powparstau, x=freqms[bb:ee],weights=1/(lmfittausstds[bb:ee]**2))
##    spec[i] = -powouttau.best_values['exponent']
#
#plt.show()
#
#if meth in 'onedim': 
#    for i in range(numfig10 +(npch/sp+1) +numFig +3):
#    ###fsuf = 'shifted'
#    ##for k in range():
#        k = numfig10 +(npch/sp+1) +numFig - i ##reverse the order
#        Summaryplot = '%s_%s_%s_%d.png'  % (pulsar,datac,meth, k+1)
#        picpathtau = newpath
#        fileoutputtau = os.path.join(picpathtau,Summaryplot)
#        plt.savefig(fileoutputtau, dpi=150)
#        print 'Saved %s in %s' %(Summaryplot,newpath)
#        plt.close()


"""Save txt files"""
    
txtpath = r'./FitTxtfiles_FluxErr_Iso'
if not os.path.exists(txtpath):
    os.makedirs(txtpath)

#1. Tau and alpha and power law amplitude -values
np.savetxt('%s/%s__%s_%s_tau_andstd.txt' %(txtpath,pulsar,datac,meth),np.vstack((taussec_highsnr,lmfitstds_highsnr*pulseperiod/nbins)))
####np.savetxt('%s/%s_%s__%s_taustd.txt' %(txtpath,pulsar,datac,meth),lmfitstds_highsnr*pulseperiod/nbins)
np.savetxt('%s/%s_%s_%s_alpha_andstd.txt' %(txtpath,pulsar,datac,meth),np.vstack((specfitdata_highsnr,specdata_err_highsnr)))
np.savetxt('%s/%s_%s_%s_powerlawamp_andstd.txt' %(txtpath,pulsar,datac,meth),np.vstack((spec_amp,spec_err_amp)))

#2. Flux-vales
np.savetxt('%s/%s_%s_%s_fluxclimbandstd.txt' %(txtpath,pulsar,datac,meth),np.vstack((fluxes_highsnr,climb_highsnr,sigmaFlux)))

#print '%s/%s_%s_%s_fluxclimbandstd.txt' %(txtpath,pulsar,datac,meth)
##3. Reduced Chi squared
np.savetxt('%s/%s_%s_%s_chi.txt' %(txtpath,pulsar,datac,meth),redchis)
##4. Fitparameters
np.savetxt('%s/%s_%s_%s_params.txt' %(txtpath,pulsar,datac,meth),bestpT_highSNR)
np.savetxt('%s/%s_%s_%s_params_std.txt' %(txtpath,pulsar,datac,meth),bestpT_std_highSNR)
##5. Freq Channels
np.savetxt('%s/%s_%s_%s_freqch.txt' %(txtpath,pulsar,datac,meth),freqms_highsnr)
#
print "Txt files saved"
#
#amp = powouttau_highsnr.best_values['amplitude']
#expp = specfitdata_highsnr
# 
#
#
#logpath = r'./LogFiles_FluxErr'
#if not os.path.exists(logpath):
#    os.makedirs(logpath)
#    
#logfilename = '%s_log_%s_%s.txt' %(pulsar,datac,meth)
#
#for k in range(11):
#    dri.log_writer(logpath,logfilename, eval('print{0}'.format(k)))
#
#print "Log files saved"

"""The next section renames arrays based on the method used, so that iso and aniso can easily be overplotted"""

"""Required arrays are"""
    #PROFILES
#
#if meth is 'onedim':
#    dataonedim = data_highsnr
#    modelonedim = model_highsnr
#    freqonedim = freqMHz_highsnr
#    tauonedim = taus_highsnr
#    tauerronedim = lmfitstds_highsnr

#plt.close('all')