# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 09:11:24 2016

@author: marisa
"""

import numpy as np


def read_header(filepath):
    f = open(filepath)
    lines = f.readlines()
    header0 = lines[0]
    h0_lines = header0.split()
    file_name = h0_lines[1]
    pulsar_name = h0_lines[3]
    nsub = int(h0_lines[5])
    nch = int(h0_lines[7])
    npol = int(h0_lines[9])
    nbins = int(h0_lines[11])
    rms = float(h0_lines[13])
#    return file_name, pulsar_name, nsub, nch, npol, nbins, rms
    return pulsar_name, nch, nbins, rms
    
    
def read_data(filepath, profilenumber, nbins):
    d = open(filepath)
    lines = d.readlines()
    
    profile_start = 2+profilenumber*(nbins+1)
    profile_end = profile_start + nbins
    
    lines_block = lines[profile_start:profile_end]
    freqc = float(lines[profile_start-1].split()[5])
    bw = float(lines[profile_start-1].split()[7])
    freqm = 10**((np.log10(freqc+ bw/2.)+ np.log10(freqc - bw/2.))/2)    
    datalist = []
    for i in range(nbins):
        data= float(lines_block[i].split()[3])
        datalist.append(data)
                
    return np.array(datalist), freqm

    

