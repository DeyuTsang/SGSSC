function [S] = SCJL1221(X,lambda)
% \lambda\|Z\|_*+\|X-XZ\|_1
%   此处提供详细说明

mu = 1;
rho = 1.01;
[d,n]=size(X);
%% Initialization
XX = X'*X;
% [U,S,V] = svd(X,'econ');
% Zt      = V(:,1:50)*V(:,1:50)';
Z = zeros(n,n);%inv(XX+lambda*eye(n,n))*XX;
S = Z;
E  = zeros(d,n);%X-X*Z;
maxIter = 200;
G1 = zeros(d,n);%;
G2 = zeros(n,n);%S-Z;
o  = ones(n,1);
eta=1e6; 
%% Iteration
for iter = 1:maxIter

    
%% Updating Z

Z = inv(mu*(XX+eye(n,n))+eta*o*o')*(mu*(XX-X'*E+S)+X'*G1-G2+eta*o*o');
%% Updating S
C  = Z+G2/mu;
S = max(sqrt( diag(C*C'))-lambda/mu,0).*C./sqrt( diag(C*C'));
%% Updating E
D  =  X-X*Z+G1/mu;
% E  =  sign(D).*max(abs(D)-1/mu,0);
E  =  max(sqrt( diag(D'*D)')-1/mu,0).*D./sqrt( diag(D'*D)');
% E  = (mu/(2*lambda+mu))*D;

%% Updating G1,G2,mu
G1 = G1+mu*(X-X*Z-E);
G2 = G2+mu*(Z-S);
mu = mu*rho;

%% check converage
err1  = norm(X-X*Z-E,'inf');
err2  = norm(Z-S,'inf');
fprintf('Iter = %.0f,err1 =%.8f,err2= %.8f \n',iter,err1,err2);


end
end
