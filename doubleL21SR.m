function [Z] = doubleL21SR(X,lambda)
% \lambda(\|Z^T\|_21+\|Z\|_21)+\|X-XZ\|_1
%   此处提供详细说明

mu = 1;
rho = 1.01;
[d,n]=size(X);
%% Initialization
XX = X'*X;
% [U,S,V] = svd(X,'econ');
% Zt      = V(:,1:400)*V(:,1:400)';
% Zt = inv(XX+lambda*eye(n,n))*XX;
Zt = zeros(n,n);%inv(XX+lambda*eye(n,n))*XX;
Z  = Zt;
S = Zt;
P = Zt;
E  = zeros(d,n);%X-X*Zt;
eta = 1e6;
maxIter = 200;
G1 = zeros(d,n);%;
G2 = zeros(n,n);
G3 = zeros(n,n);
o  = ones(n,1);
%% Iteration
for iter = 1:maxIter

%% Updating Z
A = X - E + G1/mu;
B = S + G2/mu;
C = P + G3/mu;
% Z = inv(XX+2*eye(n))*(X'*A+B+C);%无正交约束
Z = inv(mu*(XX+2*eye(n,n))+eta*o*o')*(mu*(X'*A+B+C)+eta*o*o');%有正交约束
%% Updating S
C  = Z-G2/mu;
% for i = 1:n
%     S(i,:) = max(sqrt( C(i,:)*C(i,:)')-lambda/mu,0).*C(i,:)/sqrt( C(i,:)*C(i,:)');
% end
S = max(sqrt( diag(C*C'))-lambda/mu,0).*C./sqrt( diag(C*C'));
%% Updating P
C = Z -G3/mu;
P = max(sqrt( diag(C'*C)')-lambda/mu,0).*C./sqrt( diag(C'*C)');



% for t = 1:100
% A      = X-E+G1/mu;
% B      = S+G2/mu;
% zz =  diag(Zt'*Zt).^(-0.5);
% Grad  = lambda*zz'.*Zt + mu*(XX*Zt-X'*A+Zt-B);
% Z  =  Zt - (eta/1.1^(t-1))*Grad;
% err = norm(Zt-Z,'fro');
% if mod(t,20)==0
% fprintf('Ziter=%.0f,err=%f\n',t,err)
% end
% if err<=1e-6
%     fprintf('Z函数收敛！！！\n')
%     break
% end
% Zt = Z;
% 
% % Zt = Zt-diag(diag(Zt));
% end



%% Updating E
D  =  X-X*Z+G1/mu;
E  =  sign(D).*max(abs(D)-1/mu,0);
% E  =  max(sqrt( diag(D'*D)')-1/mu,0).*D./sqrt( diag(D'*D)');
% E  = (mu/(2*lambda+mu))*D;

%% Updating G1,G2,mu
G1 = G1+mu*(X-X*Z-E);
G2 = G2+mu*(S-Z);
G3 = G3+mu*(P-Z);
mu = mu*rho;

%% check converage
err1  = norm(X-X*Z-E,'inf');
err2  = norm(Z-S,'inf');
err3  = norm(P-Z,"inf");
fprintf('Iter = %.0f,err1 =%.8f,err2= %.8f,err3= %.8f \n',iter,err1,err2,err3);


end
end
