function JVT_result=JVT(datasets,CoMa,al,norm,t)

%% Code designed and realised by Carlo Vittorio Cannistraci on 27/05/2022 
% and revised by Aga Syed Shafat Ali on 30/05/2022

%% help

%INPUT
% --> datasets: reports in a cell array the name of datasets to load, for instance:
% dataset={'carcinoma_fantom5_posHp1','developing_neurons_posHy1','MCF10A_posHy1'}; 
% in case of the tree datasets of posHp1 of the article Urbanska et al. 
% each dataset should be organized with four variables as follow:
% x: is the dataset with features (genes) or rows and samples on columns
% common_gene_names: is a cell array the gene name list that is equal for each dataset
% because all the datasets should be organized having the same genes and in
% the same order on the rows
% sample_lables: is a cell array reporting for each sample the name of one
% of the two groups to which it belongs in the dataset (this is different
% for each dataset)
% AUC_label: is a character variable with the name of one of the two groups
% of samples that should be used as positive class for AUC-ROC evaluation. 

% --> CoMa: is a cell array with the names of the genes used for computing
% the Combinatorial Marker (CoMa)for instance:
% CoMa={'cav1','fhl2','igfbp7','tagln','thbs1'}; in the article Urbanska et al. 

% --> al: is a numerical array of the same length of CoMa and that contains only 1 and -1 values
% respectively for each gene in the CoMa variable.
% 1 means that the respective gene does not need direction alignment before
% compression.
% -1 indicates that the repective gene needs direction alignment before compression
% according to the combinatorial marker formula proposed by the user. 
% For instance: al=[1,1,1,1,1];  in the article Urbanska et al. 

% --> norm: is a character variable that indicates the type of
% gene normalization which is applied over the features (genes which are on rows)
% of the datasets. This is important to avoid that, for a certain dataset, when we compress
% the CoMa genes together in a unique marker, the expression of some genes can
% dominate on others because of their magnitude. 
% The user can change the code as for its preference but here we offer two
% options: 
% 'log10': the function used is x=log10(1+x) to avoid that gene expression
% zero offers -Inf result
% 'zscore': a z-score normalization that removes the mean and divided for
% the standard deviation each fatures(genes on the rows).
% For instance: norm='log10'; in the article Urbanska et al. 
% if the inserted normalization does not match with one of the two
% available, the code uses as default norm='log10'

% --> t: is a numerical (integer) variable that indicates the number of
% sampling executed to build the null model distribution.

%OUTPUT
% --> JVT_result: is an output variable with the table with the results. I advice to call the function
% JVT.mat without any ';' at the end of the instruction because in this way
% you will be able to visualize the table directly on your 'Command Window' in MATLAB
% --> JVT_result.xls: is an output file that reports the result table. 
% It is automatically created in the same folder where this function is
% located. Pay attention that if you call this function different times with different
% input variables, the JVT_results.xls will be overwritten at each of
% these times and will report only the results of the last call. 

%% initialization

dt=length(datasets);
load(datasets{1}); cgl=length(common_gene_names); 
cm=length(CoMa); 

ict=zeros(1,cm); %location of the genes in the gene name list
for i=1:cm
   ict(i) = find(strcmp(common_gene_names,CoMa(i))); 
end


%% create null distribution

%%%%%%%%%% start initialization of variables for single gene null distribution
li=1:cgl; 
te=setdiff(li,ict); %returns the data in li that is not in ict, 
f=length(te);
li=[te ict];  
if  t < (cgl-cm) ||  t > (cgl-cm) % if t == (cgl-cm) there is not need of reshuffling 
    li=[te(randi(f,[1,t])) ict];
end
lil=length(li); % t + cm times
ROC_area=zeros(lil,dt);  ROC_area_cm=zeros(t+1,dt); 

%%%%%%%%%% finish initialization of variables for single gene null distribution

for k=1:dt % for each dataset
    
    load(datasets{k})  
    if strcmp(norm,'log10')
        x=log10(1+x);
    elseif strcmp(norm,'zscore')
        x=zscore(x,[],2);
    else 
        if k==1
    disp('########################################################################################')        
    disp('# WARNING: inserted normalization is not available. Hence, we used as default: "log10" #')  
    disp('########################################################################################')        
        end
    x=log10(1+x);    
    end
    
    
 for i=1:lil %for t + cm times

 %%% start single gene %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 [~,~,~,AUCR]=perfcurve(sample_labels,x(li(i),:),AUC_label); %instead of the AUC-ROC 
 % another performance measure can be adopted in relation with the needs of
 % the used
  if AUCR<0.5
    AUCR=1-AUCR; 
  end
 ROC_area(i,k)=AUCR;
 %%% finish single gene %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%% start combinatorial marker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i<=(t+1) 

   if i<=t
    rcm=te(randi(f,[1,cm])); %random combinatorial marker gene list indices
   elseif i==(t+1) 
    rcm=ict;
   end
   
   %%% start direction alignment
   xt=x(rcm,:); 
   for sd=1:cm 
   if al(sd)==-1
       mew=mean(xt(sd,:));
   xt(sd,:)=2*mew - xt(sd,:);   
   end   
   end
   %%% finish direction alignment
   xt=mean(xt,1);
   
  [~,~,~,AUCRcm]=perfcurve(sample_labels,xt,AUC_label); %instead of the AUC-ROC 
 % another performance measure can be adopted in relation with the needs of
 % the used
  if AUCRcm<0.5
    AUCRcm=1-AUCRcm; 
  end
 ROC_area_cm(i,k)=AUCRcm;  
 end
%%% finish combinatorial marker %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 end
end
ROC_area=min(ROC_area,[],2); % instead of the minimum the mean operator
% can be used if a more relaxed estimation is wished
ROC_area_cm=min(ROC_area_cm,[],2); % instead of the minimum the mean operator
% can be used if a more relaxed estimation is wished

%% calculate AUC-ROC and JVT P-value

result=zeros(2,cm+1);
for i=1:cm
result(1,i) = ROC_area(li==ict(i));
result(2,i) = sum(ROC_area>=result(1,i))/lil;
end
result(1,cm+1) = ROC_area_cm(t+1);
result(2,cm+1) = sum(ROC_area_cm>=result(1,cm+1))/(t+1);

% visualize the final results
CoMa{cm+1}='CombMarker';
JVT_result=array2table(result,'VariableNames',CoMa,'RowNames',{'min AUC-ROC','JVT p-Value'});
writetable(JVT_result,'JVT_result.xls','WriteRowNames',1)

