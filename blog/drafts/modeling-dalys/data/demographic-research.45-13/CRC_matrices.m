%Created by Hal Caswell and Silke van Daalen
% modified June 2021
% to accompany:
% H. Caswell and S.F. van Daalen, 2021, Healthy longevity from
% incidence-based models: More kinds of health than stars in the sky.
% Demographic Research 45:397-452
%
% Analysis of the colorectal cancer model of Wu et al. 2006, BMC Cancer 6:136
% This program creates and stores the multistate Markov chain matrices
% based on information in Wu et al. (2006)

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

%construct the transition matrix from mortality and transition rates

%ages
om=50; %age range (50 to 100)

%mortality rates
%from Human Mortality Database, Taiwan, 2002, sexes combined

load taiwan_mortality.txt;
mu=taiwan_mortality(50:50+om); %mortality starting at age 50

%transition rates
% lambda_1 = transition rate from normal cells to small adenoma (age dependent)
t=[50 55 60 65 70];
lambda_1=[0.00836  0.00990  0.01156  0.01333  0.01521]; %from Wu et al
%interpolate to get single years of age
temp=interp1(t,lambda_1,50:50+om-1,'linear','extrap');
lambda_1(1:om)=temp;

%other rates (age independent), from Wu et al
lambda=[3.46e-2
    2.15e-2
    3.697e-1
    2.382e-1
    4.852e-1
    3.02e-2
    2.099e-1];

%some useful matrices and numbers
% number of transient stages
tau=7;
% number of absorbing states
alpha=2;
%total stages
s=tau+alpha;
% total combinations ofages and stages
stot=tau*om;
% identity matrices
Itau=eye(tau);
Iom=eye(om);
Ialpha=eye(alpha);

%%%% calculation of age-dependent stage transition matrices

% create the intensity matrix Q for each age class
 for i=1:om
    Q{i}(2,1)=lambda_1(i);
    Q{i}(3,2)=lambda(1);
    Q{i}(4,3)=lambda(2);
    Q{i}(5,4)=lambda(3);
    Q{i}(6,4)=lambda(4);
    Q{i}(7,5)=lambda(5);
    Q{i}(8,6)=lambda(6);
    Q{i}(8,7)=lambda(7);
    
    Q{i}(9,1:7)=mu(i);
    Q{i}(8,8)=0;
    Q{i}(9,9)=0;
    
    Q{i}(:,:)=Q{i}(:,:)-diag(sum(Q{i}(:,:)));
    
    % create the Markov chain transition matrix using matrix exponential
    P{i}=expm(Q{i});
    % extract the transient matrix
    U{i}=P{i}(1:tau,1:tau);
    % extract the mortality matrix
    M{i}=P{i}(tau+1:tau+alpha,1:tau);
    
 end %for i

%%% create multistate matrices

% subdiagonal age transition matrix
D=diag(ones(om-1,1),-1);

%block diagonal age transition matrix
bbD=kron(Itau,D);

%block diagonal transient matrix
bbU=blkdiag(U{:});

%vec permutation matrix
K=vecperm(tau,om);

%the age-stage multistate transient matrix
Utilde=K'*bbD*K*bbU;

%age-stage multistate mortality matrix
 Mtilde=cat(2,M{:});

% age-stage multistate Markov chain transition matrix
Ptilde=[Utilde zeros(tau*om,2);
    Mtilde eye(alpha)];

%fundamental matrix of the multistate Markov chain
Ntilde=inv(eye(stot)-Utilde);

%%%% matrix construction completed

%save the matrices in a .mat file
save('CRCmatrices','Ptilde', 'Utilde', 'Mtilde', 'Ntilde', 'K', 'bbD', 'om', 'tau', 'alpha') 
