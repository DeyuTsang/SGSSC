%--------------------------------------------------------------------------
% [C,A] = lrsc_noisy(D,tau,alpha,useSampledPolynoimial)
% Low Rank Subspace Clustering algorithm for noisy data lying in a 
% union of subspaces
%
% (C,A) = argmin |C|_* + tau/2*|A - AC|_F^2 + alpha/2*|D - A|_F^2 s.t. C = C'
%
% C = affinity matrix
% A = clean data matrix whose columns are points in a union of subspaces
% D = noisy data matrix whose columns are points in a union of subspaces
% tau = scalar parameter 
% alpha = scalar parameter 
% useSampledPolynomial = when a 4th argument is passed, the code used
%                        a discretization of the polynomial thresholding
%                        function
%--------------------------------------------------------------------------
% Copyright @ Rene Vidal, November 2012
% Edited by Paolo Favaro, December 2012
%--------------------------------------------------------------------------

function [A,C] = lrsc_noisy(D,tau,alpha,~)

if nargin < 3
    alpha = 0.5*tau;
end
if nargin < 2
    tau = 100/norm(D)^2;
    alpha = 0.5*tau;
end

if nargin<=3
    [A,C] = polythreshCF(D,tau,alpha);
else
    [A,C] = polythreshAP(D,tau,alpha);
end