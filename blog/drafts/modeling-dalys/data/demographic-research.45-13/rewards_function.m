%Created by Hal Caswell and Silke van Daalen
%function to compute mean, variance, and skewness of reward
% modified June 2021
%
% % code to accompany:
% H. Caswell and S.F. van Daalen, 2021, Healthy longevity from
% incidence-based models: More kinds of health than stars in the sky.
% Demographic Research 45:397-452

function out=rewards_function(P,U, R1,R2,R3,R4)

% inputs
% P = markov chain transition matrix
% U = markov chain transient matrix
% R1, R2, R3, R4 = matrices of first four moments (about the origin) 
% of rewards associated with each transition
%
% output is a structure with moments, statistics, and other useful
% information


%sizes of matrices
%tau = number of transient states
%alpha = number of absorbing states
%s=total number of states
[s,s]=size(P);
[tau,tau]=size(U);
alpha=s-tau;

%vector of ones
e=ones(s,1);

% matrix Z to remove absorbing states
Z = [eye(tau) zeros(tau,alpha)];

% identity matrix
Itau= eye(tau);

%fundamental matrix
N=inv(Itau-U);

% vectors of moments of lifetime rewards

rho1 = N'*(Z*(((P.*R1).')*e));

rho2 = N'*(Z*(((P.*R2).')*e)+2*Z*((P.*R1).')*(Z.')*rho1);

rho3 = N'*(Z*(((P.*R3).')*e)+3*Z*((P.*R2).')*(Z.')*rho1 + ...
    3*Z*((P.*R1).')*(Z.')*rho2);

rho4 = N'*(Z*(((P.*R1).')*e)+4*Z*((P.*R3).')*(Z.')*rho1 + ...
    6*Z*((P.*R2).')*(Z.')*rho2 +4*Z*((P.*R1).')*(Z.')*rho3);

%statistics calculated from the moment vectors

%variance and standard deviation
var=rho2 - rho1.^2;
std=sqrt(var);

%coefficient of variation
%use pseudoinverse to avoid division by zero
cv=pinv(diag(rho1))*std;

%skewness
%use the pseudo-inverse
skew = (pinv(diag(var))^(3/2))*(rho3 - 3*rho2.*rho1 + 2*(rho1.^3));

%kurtosis
kurt=(-3*rho1.^4 + 6*rho2.*rho1.^2 - 4*rho1.*rho3 + rho4)./(var.^2) - 3;

% k-fold differences (lifetime)
% k=2
 arg=-log(2)./(sqrt(2*log(cv.^2+1)));
 k2 = 2*normcdf(arg);
 
% k=10
arg=-log(10)./(sqrt(2*log(cv.^2+1)));
 k10 = 2*normcdf(arg);

%create output structure returned by the function
out.rho1=rho1;
out.rho2=rho2;
out.rho3=rho3;
out.rho4=rho4;
out.var=var;
out.std=std;
out.cv=cv;
out.skew=skew;
out.kurt=kurt;
out.k2=k2;
out.k10=k10;
out.P=P;
out.R1=R1;
out.R2=R2;
out.R3=R3;




