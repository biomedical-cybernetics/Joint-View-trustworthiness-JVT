# install the necessary packages first, using:
# "pip install <package_name>" (for pip users) or 
# "conda install <package_name>" (for conda users)
import scipy.io # for loading datasets in .mat format
import numpy as np # for handling arrays, matrices and relevant computation
import pandas as pd # for result handling and saving
from sklearn.metrics import roc_auc_score # for generating AUC-ROC

def JVT(datasets,CoMa,al,norm,t):    


# Code designed by Carlo Vittorio Cannistraci on 27/05/2022
# Code tested and converted in Python by Syed Shafat Ali on 04/06/2022

# help
# INPUT
# --> datasets: reports in a cell array the name of datasets to load, for instance:
# dataset={'carcinoma','developing_neurons','MCF10A'}; in case of the three datasets of
# posHp1 of the article Urbanska et al. 
# each dataset should be organized with four variables as follow:
# x: is the dataset with features (genes) or rows and samples on columns
# common_gene_names: is a cell array the gene name list that is equal for each dataset
# because all the datasets should be organized having the same genes and in
# the same order on the rows
# datasets created in MATLAB, but do no require any special conversion before using
# in Python. This Python code handles the conversion by itself
# sample_lables: is a cell array reporting for each sample the name of one
# of the two groups to which it belongs in the dataset (this is different
# for each dataset)
# AUC_label: is a character variable with the name of one of the two groups
# of samples that should be used as positive class for AUC-ROC evaluation. Not 
# specifically required in our Python code
# --> CoMa: is a cell array with the names of the genes used for computing
# the Combinatorial Marker (CoMa)for instance:
# CoMa={'cav1','fhl2','igfbp7','tagln','thbs1'}; in the article Urbanska et al. 
# --> al: is a numerical array of the same length of CoMa and that contains only 1 and -1 values
# respectively for each gene in the CoMa variable.
# 1 means that the respective gene does not need direction alignment before
# compression.
# -1 indicates that the repective gene needs direction alignment before compression
# according to the combinatorial marker formula proposed by the user. 
# For instance: al=[1,1,1,1,1];  in the article Urbanska et al. 

# --> norm: is a character variable that indicates the type of
# gene normalization which is applied over the features (genes which are on rows)
# of the datasets. This is important to avoid that, for a certain dataset, when we compress
# the CoMa genes together in a unique marker, the expression of some genes can
# dominate on others because of their magnitude. 
# The user can change the code as for its preference but here we offer two
# options: 
# 'log10': the function used is x=log10(1+x) to avoid that gene expression
# zero offers -Inf result
# 'zscore': a z-score normalization that removes the mean and divided for
# the standard deviation each fatures(genes on the rows).
# For instance: norm='log10'; in the article Urbanska et al. 
# if the inserted normalization does not match with one of the two
# available, the code uses as default norm='log10'

# --> t: is a numerical (integer) variable that indicates the number of
# sampling executed to build the null model distribution.

# OUTPUT
# --> JVT_result: is an output variable with the table with the results. 
# --> JVT_result.xlsx: is an output file that reports the result table. 
#  It is automatically created in the same folder where this function is
#  located. Pay attention that if you call this function different times with different
#  input variables, the JVT_result.xls will be overwritten at each of
#  these times and will report only the results of the last call. 

#initialisation




    dt=len(datasets)
    ds = scipy.io.loadmat("E:/Github final/data/"+datasets[0]+'.mat')
    common_gene_names=[]
    for com_gen in ds['common_gene_names']:
        common_gene_names.append(list(com_gen[0]))
    common_gene_names = [cg for com_gen in common_gene_names for cg in com_gen]
    cgl=len(common_gene_names)
    cm=len(CoMa)
    ict=[0]*cm
    for i in range(cm):
        ict[i]=common_gene_names.index(CoMa[i])

    # create null distribution
    ## start initialization of variables for single and comb gene null distribution
    li=set(np.arange(cgl))
    te=list(li.difference(set(ict))) #returns the data in li that is not in ict
    f=len(te)
    li=te
    li.extend(ict)

    if t<(cgl-cm) or t>(cgl-cm): # if t == (cgl-cm) there is not need of reshuffling
        rand=np.round((np.random.uniform(1,f,t)))
        rand=list(map(int, rand))
        li=[te[i] for i in rand]
        li.extend(ict)
    lil=len(li)
    ROC_area=np.zeros((lil,dt))
    ROC_area_cm=np.zeros(((t+1),dt))
    ## finish initialization of variables for single and comb gene null distribution

    t=t-1 # required for Python indexing which starts from 0

    for k in range(dt): # for each dataset
        ds = scipy.io.loadmat("E:/Github final/data/"+datasets[k]+'.mat')
        x=ds['x'] # the actual data
        sample_labels=[]
        r=ds['sample_labels'].shape[0]
        c=ds['sample_labels'].shape[1]
        if r==1:
            ds['sample_labels']=ds['sample_labels'].reshape(c,r)
        for SL in ds['sample_labels']:
            sample_labels.append(list(SL[0]))
        sample_labels = [sample_label for SL in sample_labels for sample_label in SL]
        AUC_label=ds['AUC_label'] # not required in Python
        if norm=='log10':
            x=np.log10(1+x)
        elif norm=='zscore':
            x=scipy.stats.zscore(x)
        else:
            if k==1:
                print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
                print('% WARNING: inserted normalization is not available. Hence, we used as default: log10 %')
                print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            x=np.log10(1+x)
        
        for i in range(lil): # for t+cm times

            ## start single gene ##
            AUCR=roc_auc_score(sample_labels,x[li[i],:])
            # instead of the AUC-ROC another performance measure can be adopted 
            # in relation with the needs of the user
            ROC_area[i,k]=max(AUCR,1-AUCR)
            #### finish single gene ####

            #### start combinatorial marker ####
            if i<=(t+1):
                if i<=t:
                    rand=np.round((np.random.uniform(1,f,cm)))
                    rand=list(map(int, rand))
                    rcm=[te[i] for i in rand]
                elif i==(t+1):
                    rcm=ict

                #### start direction alignment
                xt=x[rcm,:]
                for sd in range(cm):
                    if al[sd]==-1:
                        mew=np.mean(xt[sd,:])
                        xt[sd,]=2*mew-xt[sd,:]
                #### finish direction alignment

                xt=np.mean(xt,axis=0)
                AUCRcm=roc_auc_score(sample_labels,xt) #instead of the AUC-ROC
                # another performance measure can be adopted in relation with 
                # the needs of the user

                ROC_area_cm[i,k]=max(AUCRcm,1-AUCRcm)
                #### finish combinatorial marker
    ROC_area=np.min(ROC_area,axis=1) # instead of the minimum, the mean operator
    # can be used if a more relaxed estimation is wished           
    ROC_area_cm=np.min(ROC_area_cm,axis=1) # instead of the minimum, the mean operator
    # can be used if a more relaxed estimation is wished


    #### calculate AUC-ROC and JVT P-value

    result=np.zeros((2,(cm+1)))
    for i in range(cm):
        result[0,i]=ROC_area[li.index(ict[i])] 
        result[1,i] = sum(ROC_area>=result[0,i])/lil
    result[0,cm] = ROC_area_cm[t+1]
    result[1,cm] = sum(ROC_area_cm>=result[0,cm])/(t+2)
    # t+2, restored the earlier deduction in t

    #### visualize the final statistic results of mean and minimum
    CoMa.append('Comb Marker')
    JVT_result=pd.DataFrame(result)
    JVT_result.columns=CoMa
    JVT_result.index=['min AUC-ROC','JVT p-value']
    JVT_result.to_excel("JVT_result.xlsx")  
    print(JVT_result)
