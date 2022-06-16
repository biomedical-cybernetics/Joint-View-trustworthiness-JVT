# set your working directory the R Code folder
setwd("E:\\Github final\\R\\")

#initialization of the input variables for positive hypothesis 2 in Table 4 of Urbanska et al
datasets = c('carcinoma_CCLE_micro_posHy2','carcinoma_CCLE_RNASeq_posHy2','carcinoma_Klijn_posHy2')
CoMa=c('cav1','fhl2','igfbp7','tagln','thbs1')
al=c(1,1,1,1,1)
norm='log10' # other option available is 'zscore'
t=1000 # in the study Urbanska et al. we used t=10000, but here for speed up
# the computation we use t=1000 since the results are comparable

# JVT computation for positive hypothesis 2 reproduce results in Table 4 of Urbanska et al.
source("JVT.R")
JVT(datasets,CoMa,al,norm,t)
