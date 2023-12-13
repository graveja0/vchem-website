%Created by Hal Caswell and Silke van Daalen
% modified June 2021
%
%program to calculate the 
% a) the number of transitions into clinical colorectal cancer (CRC)stages
% b) the statistics of the years of life lost to colorectal cancer
% based on the model of Wu et al. 2006, BMC Cancer 6:136
%
% code to accompany:
% H. Caswell and S.F. van Daalen, 2021, Healthy longevity from
% incidence-based models: More kinds of health than stars in the sky.
% Demographic Research 45:397-452

% stages
%1 = normal
% 2 = small adenoma
% 3 = large adenoma
% 4 = preclinical early CRC
% 5 = preclinical late CRC
% 6 = clinical early CRC
% 7 = clinical late CRC
% 8 = CRC death
% 9 = other causes of death

%ages
om=50; %age range (50 to 100)

%load file with matrices
load CRCmatrices


%NUMBER OF TRANSITIONS INTO CLINICAL CANCER
% Healthy longevity = count of transitions into stages 6 and 7
% over entire lifetime

B=zeros(tau,tau,om);
B(6,4,:)=1;
B(7,5,:)=1;

C=zeros(alpha,tau,om);

ages=zeros(om,1);
ages(1:om)=1;

%create block diagonal transition matrix
bbB=kron(diag(ages),B(:,:,1));
Btilde=K'*bbD*K*bbB;
%create block mortality matrix
Ctilde=kron(ones(1,om),C(:,:,1));

%create reward matrices
R1=[Btilde zeros(tau*om,alpha);
    Ctilde zeros(alpha,alpha)];
% the number of transitions is a fixed reward
R2=R1.*R1;
R3=R1.*R2;
R4=R1.*R3;

%call rewards function
out_allages=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);

%YEARS OF LIFE LOST DUE TO COLORECTAL CANCER
%calculating life lost due to transitions into death due to CRC

%moments of remaining longevity
eta1=(sum(Ntilde))';
eta2=(eta1'*(2*Ntilde - eye(tau*om)))';
eta3=(eta1'*(6*Ntilde^2 -6*Ntilde +eye(tau*om)))';
eta4=(eta1'*(24*Ntilde^3 + 36*Ntilde^2 + 14*Ntilde - eye(tau*om)))';

%reward matrices based on transition to stage 8 (CRC deatsh)

Btilde=zeros(tau*om);
Ctilde=zeros(alpha,tau*om);
Ctilde(1,:)=eta1';
R1=[Btilde zeros(tau*om,alpha);
    Ctilde zeros(alpha,alpha)];

Ctilde(1,:)=eta2';
R2=[Btilde zeros(tau*om,alpha);
    Ctilde zeros(alpha,alpha)];


Ctilde(1,:)=eta3';
R3=[Btilde zeros(tau*om,alpha);
    Ctilde zeros(alpha,alpha)];

Ctilde(1,:)=eta4';
R4=[Btilde zeros(tau*om,alpha);
    Ctilde zeros(alpha,alpha)];

%call rewards function
out_YLL=rewards_function_sept2020(Ptilde,Utilde,R1,R2,R3,R4);




