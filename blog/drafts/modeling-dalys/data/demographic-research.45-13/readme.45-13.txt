This folder contains files that provide the Matlab code for calculations of the examples presented in Caswell and van Daalen (Demographic Research 45:397-452, 2021) for the colorectal cancer (CRC) model of Wu et al. (2016). Matlab programs have been used under version R2020a. 


CRC_matrices.m --- a script that computes the multistate Markov chain matrices for the Wu et al. model. 

CRC_occupancy.m --- a script that calculates the healthy longevity in terms of occupancy, for the definitions used in the paper as examples. Easily modified to any desired occupancy.

CRC_transitions.m --- a script that calculates healthy longevity in terms of transitions, for the examples used in the paper. Easily modified to any desired transitions. 

reward_function.m --- a function that returns the first four moments of the rewards, however those may be defined.

vec.m -- a simple function to return the vec of an array

vecperm.m -- a function that calculates the vec-permutation matrix for specified dimensions

CRCmatrices.mat --- a Matlab data file containing all the matrices used in the analysis. This makes it unnecessary to create the matrices in the program CRC_matrices.m

taiwan_mortality.txt -- a text file that is loaded by the longevity programs to provide the age-specific mortality rates for Taiwan. The values are mortality probabilities (q_x) for age classes 1 to 110, for the year 2000, both sexes combined. Data obtained from the Human Mortality Database, https://mortality.org.




