% vecperm
% function to calculate the vec permutation matrix K of index m,n
% let X be a m x n matrix, and X' the transpose of X
% then K satisfies 
% vec(X') = K*vec(X)
% to accompany:
% H. Caswell and S.F. van Daalen, 2020, Healthy longevity from
% incidence-based models: More kinds of health than stars in the sky.
% Demographic Research 45:397-452

function K = vecperm(m,n)

K = zeros(m*n);
a = zeros(m,n);
for i = 1:m
    for j = 1:n
        e = a;
        e(i,j) = 1;
        K = K + kron(e,e');
    end
end


        
