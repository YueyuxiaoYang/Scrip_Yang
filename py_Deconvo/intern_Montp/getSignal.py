# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 21:32:19 2018

@author: kevin1024
"""
import numpy as np

def getSignal(PolyPosition, Parameters):
    [FreqEchS, FreqEchImg, DureeAnalysee, NSondeFluo, NSondeParIntensite,\
        TaillePreMarq,TailleSeqMarq, TaillePostMarq, VitessePolymerase,frame_num] = Parameters
    Intensity_for_1_Polym = NSondeFluo/NSondeParIntensite;
    T1 = round((TaillePreMarq/VitessePolymerase)*FreqEchS);
    T2 = round((TailleSeqMarq/VitessePolymerase)*FreqEchS);
    T3 = round((TaillePostMarq/VitessePolymerase)*FreqEchS);
    PolyPosition = round(PolyPosition);
    first = np.zeros(int(PolyPosition+T1));
    second = np.arange(Intensity_for_1_Polym/T2,Intensity_for_1_Polym,Intensity_for_1_Polym/T2)
    third = np.repeat(Intensity_for_1_Polym,T3);
    fourth = np.concatenate((first, second, third),axis=0);  # the assembled signal
    len_1_sig = T1+T2+T3; # how many "interval" for 1 signal
    fifth = np.array([0]);
    if(len(fourth)<len_1_sig+round(frame_num*(FreqEchS/FreqEchImg))+1):
        fifth = np.zeros(int((frame_num*(FreqEchS/FreqEchImg))+1));  # the assembled signal is smaller than the desired simulation duration add some 0 to the end
    sixth = np.concatenate((fourth,fifth),axis=0);
    #     sixth = sixth(1:round(frame_num*FreqEchS/FreqEchImg));  # In all cases keep only the signal of the simulated duration length
    #     signal = sixth(len_1_sig:round(FreqEchS/FreqEchImg):len_1_sig+frame_num*round(FreqEchS/FreqEchImg));  # return only the values corresponding to the experimental value times points
    ii=1;
    n_signal = np.zeros(int(frame_num));
    while ii <= frame_num:
        index = int((ii-1)*(FreqEchS/FreqEchImg)+len_1_sig); # image start at len_1_sig, because 
        n_signal[ii-1] = sixth[index-1];   
        ii= ii+1;         
    return n_signal

def sumSignal(Trans_positions,Parameters):
    l_signal = len(getSignal(1, Parameters));
    Sum_signals = np.zeros(l_signal);
    for posi_i in Trans_positions:
        signal = getSignal(posi_i, Parameters);
        Sum_signals = Sum_signals + signal; 
    return Sum_signals