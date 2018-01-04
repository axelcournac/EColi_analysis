# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 14:38:25 2015
@author: Axel KournaK
To smooth genomic signal (binned with restriction frags or kb) taking into account circularity of genome
"""

import numpy as np

def sliding_func(signal,  nb_frags):
    smoothed = np.zeros((len(signal), 1));
    for i in range(0,len(signal) ) :
        nf=0;
        for k in range(i-nb_frags/2,i+nb_frags/2 +1):
            kk=k;
            if kk >= len(signal) :
                kk = kk - len(signal);
            if kk <0:
                kk= len(signal) +k;
            if signal[kk] > 0:
                smoothed[i]=smoothed[i]+signal[kk];
                nf+=1;
        smoothed[i] = smoothed[i] / nf;
    
    return smoothed