%Created by Hal Caswell and Silke van Daalen
% modified June 2021
% program to calculate healthy longevity as occupancy
% of sets of health stages at specified ages
%
% code to accompany:
% H. Caswell and S.F. van Daalen, 2021, Healthy longevity from
% incidence-based models: More kinds of health than stars in the sky.
% Demographic Research 45:397-452
%
%analysis of the colorectal cancer model of Wu et al. 2006, BMC Cancer
% 6:136
%
% stages in the model
% 1 = normal cells
% 2 = small adenoma
% 3 = large adenoma
% 4 = preclinical early CRC
% 5 = preclinical late CRC
% 6 = clinical early CRC
% 7 = clinical late CRC
% 8 = death from CRC
% 9 = death from other causes

 %load the matrices for the CRC model
load('CRCmatrices');

%some useful matrices and numbers
%transient stages
tau=7;
%absorbing states
alpha=2;
%total stages
s=tau+alpha;
% total combinations
stot=tau*om;
%matrices
Itau=eye(tau);
Iom=eye(om);
Ialpha=eye(alpha);
eentau=ones(tau,1);  %"een" is "one" in Dutch
eenom=ones(om,1);
eenalpha=ones(alpha,1);
eentot=ones(tau*om,1);

% Healthy longevity calculations. Several different definitions are
% calculated

% H= health matrix, dimensions tau X om (i.e., stage x age)
% note: in Matlab, ~X is the logical complement of X

% TOTAL LONGEVITY
% Healthy longevity = time spent in any living stage  

H=ones(tau,om);
h=vec(H);

B1= h*h' +0.5*h*(~h') +0.5*(~h)*h';
C1=0.5*eenalpha*h';

R1=[B1 zeros(tau*om,alpha);
    C1 zeros(alpha,alpha)];

%occupancy is a fixed reward, so...
R2=R1.^2;
R3=R1.^3;
R4=R1.^4;

%call the rewards function
out_long=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);


%NORMAL CELL LONGEVITY
%Healthy longevity = time spent in the normal cell stage (stage 1)

H=zeros(tau,om);
H(1,:)=1;
h=vec(H);
 
B1= h*h' +0.5*h*(~h') +0.5*(~h)*h';
C1=0.5*eenalpha*h';

R1=[B1 zeros(tau*om,alpha);
    C1 zeros(alpha,alpha)];

%occupancy is a fixed reward, so
R2=R1.^2;
R3=R1.^3;
R4=R1.^4;

%call the rewards function
out_normal=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);


% EARLY AGE CANCER
% Healthy longevity = occupancy of any stage with cancer (4-7), at ages <= 65

H=zeros(tau,om);
H([4 5 6 7],1:15)=1;
h=vec(H);

B1= h*h' +0.5*h*(~h') +0.5*(~h)*h';
C1=0.5*eenalpha*h';

R1=[B1 zeros(tau*om,alpha);
    C1 zeros(alpha,alpha)];

%occupancy as a fixed reward, so
R2=R1.^2;
R3=R1.^3;
R4=R1.^4;

%call the rewards function
out_CRCearly=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);

% LATE AGE CANCER
% Healthy longevity = occupancy of any of the stages (4-7) with cancer after age 65

H=zeros(tau,om);
H([4 5 6 7],16:end)=1;
h=vec(H);

B1= h*h' +0.5*h*(~h') +0.5*(~h)*h';
C1=0.5*eenalpha*h';

R1=[B1 zeros(tau*om,alpha);
    C1 zeros(alpha,alpha)];

%occupancy as a fixed reward
R2=R1.^2;
R3=R1.^3;
R4=R1.^4;

%call rewards function
out_CRClate=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);


% CANCER-FREE LONGEVITY
% Healthy longevity = occupancy of any of the cancer-free stages (1-3)

H=zeros(tau,om);
H([1 2 3],:)=1;
h=vec(H);

B1= h*h' +0.5*h*(~h') +0.5*(~h)*h';
C1= 0.5*eenalpha*h';

R1=[B1 zeros(tau*om,alpha);
    C1 zeros(alpha,alpha)];

%occupancy is a fixed reward
R2=R1.^2;
R3=R1.^3;
R4=R1.^4;

%call rewards function
out_cancerfree=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);



