#install the necessary packages first, using "install.package('package_name')"
library(pROC) # library for generating AUC-ROC
library(R.matlab) # library for loading datasets in .mat format
library(xlsx) # for writing the result to .xlsx file

JVT <- function(datasets,CoMa,al,norm,t){


# Code designed by Carlo Vittorio Cannistraci on 27/05/2022
# Code tested and converted in R by Syed Shafat Ali on 30/05/2022

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
# in R. This R code handles the conversion by itself
# sample_lables: is a cell array reporting for each sample the name of one
# of the two groups to which it belongs in the dataset (this is different
# for each dataset)
# AUC_label: is a character variable with the name of one of the two groups
# of samples that should be used as positive class for AUC-ROC evaluation. Not 
# specifically required in our R code
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

dt=length(datasets)
ds=readMat(paste("E:\\Github final\\data\\",datasets[1],".mat",sep=''))
common_gene_names=as.character(unlist(ds[2])) # ds[2] is the common gene names
cgl=length(common_gene_names)
cm=length(CoMa)
ict=rep(0,cm)
for(i in 1:cm){ #location of the genes in the gene name list
   ict[i] = which(common_gene_names==CoMa[i]) 
}

#create null distribution
##### start initialization of variables for single and comb gene null distribution
li=1:cgl

te=setdiff(li,ict) #returns the data in li that is not in ict, 
f=length(te)
li=append(te,ict)
if  (t < (cgl-cm) ||  t > (cgl-cm)){ # if t == (cgl-cm) there is not need of reshuffling 
    li=append(te[round(runif(t, 1 , f))],ict)
}
lil=length(li); # t + cm times
ROC_area=matrix(0,lil,dt) 
ROC_area_cm=matrix(0,t+1,dt) 

##### finish initialization of variables for single and comb gene null distribution

for(k in 1:dt){ # for each dataset
	ds=readMat(paste("E:\\Github final\\data\\",datasets[k],".mat",sep=''))
	x=ds[4]$x # ds[4]$x is actually the data x
	sample_labels=as.character(unlist(ds[3])) # ds[3] is sample labels
	AUC_label=as.character(unlist(ds[1])) # ds[1] is AUC label; not required in R
	if(norm=='log10'){
		x=log(1+x,10)
	}
	else if(norm=='zscore'){
		x = (x - mean(x)) / sd(x) # zscore normalisation
	}
	else{
			if(k==1){
				print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
				print('% WARNING: inserted normalization is not available. Hence, we used as default: log10 %')
				print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
			}
		x=log(1+x,10)
	}
	for(i in 1:lil){ # for t+cm times

		## start single gene ##
		AUCR=as.numeric(auc(roc(sample_labels,x[li[i],],quiet=T)))
		#instead of the AUC-ROC another performance measure can be adopted 
		#in relation with the needs of the user
		ROC_area[i,k]=max(AUCR,1-AUCR)
		#### finish single gene ####

		#### start combinatorial marker ####
		if(i<=(t+1)){
			if(i<=t){
				rcm=te[round(runif(cm, 1 , f))]
			}
			else if(i==t+1){
				rcm=ict
			}

			#### start direction alignment
			xt=x[rcm,]
			for(sd in 1:cm){
				if(al[sd]==-1){
				mew=mean(xt[sd,])
				xt[sd,]=2*mew-xt[sd,]
				}
			}
			#### finish direction alignment

			xt=colMeans(xt)
			AUCRcm=as.numeric(auc(roc(sample_labels,xt,quiet=T))) #instead of the AUC-ROC 
 			# another performance measure can be adopted in relation with the needs of
 			# the user

			ROC_area_cm[i,k]=max(AUCRcm,1-AUCRcm)
			#### finish combinatorial marker
		}
	}

}
ROC_area=apply(ROC_area, 1, FUN = min) # instead of the minimum, the mean operator
# can be used if a more relaxed estimation is wished
ROC_area_cm=apply(ROC_area_cm, 1, FUN = min) # instead of the minimum, the 
# mean operator can be used if a more relaxed estimation is wished


#### calculate AUC-ROC and JVT P-value

result=matrix(0,2,cm+1)
for(i in 1:cm){
	result[1,i] = ROC_area[li==ict[i]]
	result[2,i] = sum(round(ROC_area,4)>=round(result[1,i],4))/lil
}
result[1,cm+1] = ROC_area_cm[t+1]
result[2,cm+1] = sum(round(ROC_area_cm,4)>=round(result[1,cm+1],4))/(t+1)

#### visualize the final statistic results of mean and minimum
CoMa=append(CoMa,'Comb Marker')
JVT_result=as.data.frame(result)
colnames(JVT_result)=CoMa
rownames(JVT_result)=c('min AUC-ROC','JVT p-value')
write.xlsx(JVT_result,'JVT_result.xlsx',col.names = TRUE,row.names = TRUE)
print(JVT_result)
}
