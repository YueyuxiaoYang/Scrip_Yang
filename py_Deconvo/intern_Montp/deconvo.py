# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 21:45:50 2018

@author: kevin1024
"""
import numpy as np 
import os
from getSignal import *
from GA import GA

'''
Import files and read experiment data
can only read csv from 1 filePath
'''
DataFilePath ='/home/kevin1024/Desktop/intern_Montp/examples_cvs/'
file_name_list = os.listdir(DataFilePath)
DataFileName = DataFilePath+file_name_list[1]
DataExp = np.genfromtxt(DataFileName,delimiter=';',skip_header=1)
DataExpSmooth = np.transpose(DataExp)[1]  

'''
 Model and Parameters 
'''
#--------transcription--------
NbrSondeFluo = 128.0; # 128 sondes per mRNA (ref publi)
TaillePreMarq = 700.0; # 700 bases (ref publi)
TailleSeqMarq = 5800.0; # 2900 bases (ref publi)
TaillePostMarq = 1600.0; # 1600 bases (ref publi)
DureeProcessing = 100.0; #  en secondes (100s temps moyen ref publi)
# --------polymerases------- 
EspaceInterPolyMin = 30.0; # en base (40 bases)
Polym_speed = 67.0; # average speed bases par seconde (Ref publi)
ProbeByIntensitie_nb = NbrSondeFluo; # how many probes are needed to observe 1 intensity
# -------sample freqence------
NombreRandomIter = 1000; #at least 1000 but 10000 is fine for 1200 secondes et FreqEchSimu1
FreqEchImg = (1.0/3); # 1/3 image par seconde
# FreqEchSimu = 2; # nb de calcul par seconde ## Attention cette donn?e doit ?tre enti?re. -> int ; no float ! ET PROP a la dur?e des calculs
# DureeSimu = 3000; # 1200 secondes Attention cette valeur ne doit pas ?tre inf?rieure a la dur?e exp?rimentale
#CoeffNormalisation = 1
Gap4TrainDef = 2;
DureeSimu = DataExp[-1][0]+DataExp[1][0] # (s)
frame_num = DataExp.shape[0];
DureeSignal = round((TaillePreMarq+TailleSeqMarq + TaillePostMarq) / Polym_speed); # (s)
DureeAnalysee = DureeSignal + DureeSimu - round(TaillePreMarq / Polym_speed); # (s)
FreqEchSimu = 1.0/(EspaceInterPolyMin/Polym_speed); # how many interval(possible poly start position) in 1s
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed));
Parameters = [FreqEchSimu, FreqEchImg,DureeAnalysee,NbrSondeFluo,\
            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num];


'''
Generate artificial data

'''

Nbr_poly= 150    
Poly_position_art = np.random.choice(range(2900),size=Nbr_poly,replace=False)
sum_signal_art = sumSignal(Poly_position_art,Parameters)
pop_GA,stats, hof=GA(sum_signal_art,Parameters,IND_SIZE = num_possible_poly,Population_Nbr=20,Max_gen=10,Nbr_poly_estimate=170)
    
'''
pop_GA.sort(key=lambda ind:ind.fitness.values)
plot(sumSignal(np.where(pop_GA[0]==1)[0],Parameters))
plot(sum_signal_art,'red')
    
'''
    
    
    
    
    
    
    
    
    