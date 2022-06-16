%% initialization of the input variables for negative Hypothesis in Table 4 of Urbanska et al
datasets = {'carcinoma_CCLE_micro_1LuS_NegHy','carcinoma_CCLE_micro_2InS_NegHy','carcinoma_Klijn_NegHy'};
CoMa={'cav1','fhl2','igfbp7','tagln','thbs1'};
al=[1,1,1,1,1]; 
norm= 'log10'; %'zscore';
t=1000; % in the study Urbanska et al. we used t=10000, but here for speed up
%the computation we use t=1000 since the results are comaprable

%% JVT computation for negative hypothesis reproduces results in Table 4 of Urbanska et al.
tic
JVT(datasets,CoMa,al,norm,t)
toc
