# set your working directory the R Code folder
setwd("E:\\Github final\\R\\")

#initialization of the input variables for negative hypothesis in Table 4 of Urbanska et al
datasets = c('carcinoma_CCLE_micro_1LuS_NegHy','carcinoma_CCLE_micro_2InS_NegHy','carcinoma_Klijn_NegHy')
CoMa=c('cav1','fhl2','igfbp7','tagln','thbs1')
al=c(1,1,1,1,1)
norm='log10' # other option available is 'zscore'
t=1000 # in the study Urbanska et al. we used t=10000, but here for speed up
# the computation we use t=1000 since the results are comparable

# JVT computation for negative hypothesis reproduce results in Table 4 of Urbanska et al.
source("JVT.R")
JVT(datasets,CoMa,al,norm,t)
