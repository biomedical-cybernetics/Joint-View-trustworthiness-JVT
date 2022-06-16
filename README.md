# Joint-View-trustworthiness-JVT

# Definition of jointly multiview discrimination
Jointly multiview discriminative means that the minimum performance of a predictor model or a marker (single or combinatorial) on different datasets is considered as a reference of its general discriminative performance. Indeed, if the minimum performance is high then it means the predictor is general in ability to discriminate the investigated phenotype. Meanwhile, this general ability to discriminate should be specific of the investigated phenotype. If a marker displays both generality and specificity to a certain phenotype then we can expect that it is a universal and specific signature associated to that phenotype. 

# Joint-wiew trustworthiness (JVT)
JVT is a resample technique that creates a null model distribution according to which an empirical p-value is computed to evaluate the probability to sample at random a single or combinatorial marker that offers a joint multview discrimination that is equal or better than the tested marker. A low JVT p-value (0.05 significant level) means that it is rare to generate at random a joint multiview marker with performance equal or better than the tested one. 
Methodological details on JVT are provided in the Methods of the article: De novo identification of universal cell mechanics gene signatures by Urbanska et al. 
In the framework of Urbanska et al. study: the theoretical concept, algorithm and applications for JVT were invented and written by Carlo Vittorio Cannistraci with feedback by Yan Ge, Marta Urbanska and Syed Shafat Ali. 

# What is containing this repository
Data, MATLAB, R and Python code files in respective dedicated folders, specifically:
+ Data: the 9 datasets that you need to run either of the MATLAB, R or Python scripts.
+ MATLAB (written by Carlo Vittorio Cannistraci and tested by Syed Shafat Ali): JVT code and three scripts to replicate the results of positive hypothesis 1, positive hypothesis 2 and negative hypothesis in Table 4 of Urbanska et al. At the beginning of the code, there is a help section that explains how to use it.
+ Python (written by Syed Shafat Ali and tested by Yan Ge): analogous functions of the MATLAB folder. You will find both the Python scripts (.py) and Jupyter Notebook scripts (.pynb) for each of positive hypothesis 1, positive hypothesis 2 and negative hypothesis.
+ R  (written by Syed Shafat Ali and tested by Yan Ge): analogous functions of the MATLAB folder.

Note: Ensure you set the paths correctly before running the scripts. The functions were tested respectively in: MATLAB 2018a or youger, Python 3.9.4, R 4.0.3. 

# Contacts
For any problem, please contact:
Carlo Vittorio Cannistraci kalokagathos.agon@gmail.com


