function [Z,F] = LapdoubleL21SR(X,lam1,lam2,c)
% \lambda(\|Z^T\|_21+\|Z\|_21)+\|X-XZ\|_1
%   此处提供详细说明
% c  = c;
mu = 0.6;
rho = 1.01;
[d,n]=size(X);
%% Initialization
XX = X'*X;
% [U,S,V] = svd(X,'econ');
% Zt      = V(:,1:400)*V(:,1:400)';
Zt = inv(XX+eye(n,n))*XX;
Z  = Zt;
S  = Zt;
J  = Zt;
E  = X-X*Zt;
eta = 0.0001;
maxIter = 100;
G1 = X-X*Zt-E;
G2 = S-Zt;
G3 = J-Zt;
W  = abs(Z)+abs(Z');
D  = diag(diag(W));
L  = D-W;
[eV,eD] = eig(L);
F  = eV(:,1:c);
for i = 1:n
P(:,i) = sum((F(i,:)-F).^2,2);
end
%% Iteration
for iter = 1:maxIter


%% Updating S
A  = Z+G2/mu;
% for i = 1:n
%     S(i,:) = max(sqrt( C(i,:)*C(i,:)')-lambda/mu,0).*C(i,:)/sqrt( C(i,:)*C(i,:)');
% end
S = max(sqrt( diag(A*A'))-1/mu,0).*A./sqrt( diag(A*A'));
%% Updating Z
for t = 1:100
Q      = X-E+G1/mu;
G      = S-G2/mu;
H      = J-G3/mu;
% for j = 1:n
% zt     = Zt(:,j);
% zz     = zt'*zt;
% grad   = lambda*(zz)^(-0.5)*zt+mu*(XX*zt-X'*A(:,j)+zt-B(:,j));
% Z(:,j) = zt-(eta/1.1^(t-1))*grad;
% end
zz =  diag(Zt'*Zt).^(-0.5);
Grad  = zz'.*Zt + mu*(XX*Zt-X'*Q+Zt-G+Zt-H);
Z  =  Zt - (eta/1.1^(t-1))*Grad;
err = norm(Zt-Z,'fro');
if mod(t,20)==0
fprintf('Ziter=%.0f,err=%f\n',t,err)
end
if err<=1e-6
    fprintf('Z函数收敛！！！\n')
    break
end
Zt = Z;

% Zt = Zt-diag(diag(Zt));
end
%% Updating J
B      = Z+G3/mu;
J      = sign(B).*(abs(B)-lam2*P/mu);
%% Updating P
W      = abs(Z)+abs(Z');
D      = diag(diag(W));
L      = D-W;
[eV,eD] = eig(L);
F  = eV(:,1:c);
for i = 1:n
P(:,i) = sum((F(i,:)-F).^2,2);
end
%% Updating E
C  =  X-X*Z+G1/mu;
E  =  sign(C).*max(abs(C)-lam1/mu,0);

%% Updating G1,G2,mu
G1 = G1+mu*(X-X*Z-E);
G2 = G2+mu*(Z-S);
G2 = G2+mu*(Z-J);
mu = mu*rho;

%% check converage
err1  = norm(X-X*Z-E,'inf');
err2  = norm(Z-S,'inf');
err3  = norm(Z-J,'inf');
fprintf('Iter = %.0f,err1 =%.8f,err2= %.8f,err3= %.8f \n',iter,err1,err2,err3);


end
end
