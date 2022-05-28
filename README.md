# Joint-View-trustworthiness-JVT

# Definition of jointly multiview discrimination
Jointly multiview discriminative means that the minimum performance of a predictor model or a marker (single or combinatorial) on different datasets is considered as a reference of its general discriminative performance. Indeed, if the minimum performance is high then it means the predictor is general in ability to discriminate the investigated phenotype. Meanwhile, this general ability to discriminate should be specific of the investigated phenotype. If a marker displays both generality and specificity to a certain phenotype then we can expect that it is a universal and specific signature associated to that phenotype. 

# Joint-wiew trustworthiness (JVT)
JVT is a resample technique that creates a null model distribution according to which an empirical p-value is computed to evaluate the probability to sample at random a single or combinatorial marker that offers a joint multview discrimination that is equal or better than the tested marker. A low JVT p-value (<0.05 significant level) means that it is rare to generate at random a joint multiview marker with performance equal or better than the tested one. 
Methodological details on JVT are provided in the Methods of the article: De novo identification of universal cell mechanics gene signatures by Urbanska et al. 
In the framework of Urbanska et al. study: the theoretical concept, algorithm, code and applications for JVT were invented and written by Carlo Vittorio Cannistraci with feedback by Yan Ge, Marta Urbanska and Syed Shafat Ali. 

# What is containing this repository
MATLAB code files and datasets, specifically:
+ JVT code. At the beginning of the code there is an help section that explains how to use it. 
+ three scritps that you can run to replicate the results of positive hytothesis 1, positive hypothesis 2 and negative hypothesis in Table 4 of Urbanska et al. 
+ the 11 datasets that you need to run the scripts


