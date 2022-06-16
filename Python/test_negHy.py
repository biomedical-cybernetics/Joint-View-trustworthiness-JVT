#!/usr/bin/env python
# coding: utf-8

import os
os.chdir("E:/Github final/Python") 
# set your working directory the 'Python Code' folder
from importlib import reload
import JVT # import the file carrying JVT function
JVT=reload(JVT)

#initialization of the input variables for negative hypothesis in Table 4 of Urbanska et al
datasets = ['carcinoma_CCLE_micro_1LuS_NegHy','carcinoma_CCLE_micro_2InS_NegHy','carcinoma_Klijn_NegHy']
CoMa=['cav1','fhl2','igfbp7','tagln','thbs1']
al=[1,1,1,1,1]
norm='log10' # available options 'log10' and 'zscore'
t=1000 # in the study Urbanska et al. we used t=10000, but here for speed up
# the computation we use t=1000 since the results are comparable

#JVT computation for negative hypothesis reproduce results in Table 4 of Urbanska et al.
JVT.JVT(datasets,CoMa,al,norm,t)
