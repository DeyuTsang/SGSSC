% min_{A>=0, A*1=1, F'*F=I}  trace(D'*A) + r*||A||^2 + 2*lambda*trace(F'*L*F)
% written by Feiping Nie on 2/9/2014
function [y, A, evs] = CAN(X, c, k, r, islocal)
% X: dim*num data matrix, each column is a data point
% c: number of clusters
% k: number of neighbors to determine the initial graph, and the parameter r if r<=0
% r: paremeter, which could be set to a large enough value. If r<0, then it is determined by algorithm with k
% islocal: 
%           1: only update the similarities of the k neighbor pairs, faster
%           0: update all the similarities
% y: num*1 cluster indicator vector
% A: num*num learned symmetric similarity matrix
% evs: eigenvalues of learned graph Laplacian in the iterations

% For more details, please see:
% Feiping Nie, Xiaoqian Wang, Heng Huang. 
% Clustering and Projected Clustering with Adaptive Neighbors.  
% The 20th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD), New York, USA, 2014.



NITER = 30;
num = size(X,2);
if nargin < 5
    islocal = 1;
end;
if nargin < 4
    r = -1;
end;
if nargin < 3
    k = 15;
end;

distX = L2_distance_1(X,X);
%distX = sqrt(distX);
[distX1, idx] = sort(distX,2);
A = zeros(num);
rr = zeros(num,1);

% 初始化操作 A是一个类似k近邻图的东西。用的似乎是CLR这篇文章的init Graph A 的算法。
for i = 1:num
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;

% 画出初始K近邻图
figure('name','init graph Matix A '); 
imshow(A,[]); colormap jet; colorbar;
XX = X';
figure('name',' init graph ');
for i=1:200
    for j=1:200
        if(A(i,j)~=0)
            plot([XX(i,1),XX(j,1)],[XX(i,2),XX(j,2)],'--gs','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','r'); hold on;
        end
    end
end;


if r <= 0
    r = mean(rr);
end;
lambda = mean(rr);      % initialize lambda by r

% 构造拉普拉斯矩阵，
% 假设A固定，求解F
A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, temp, evs]=eig1(L0, c, 0);%得到两百个特征向量

if sum(evs(1:c+1)) < 0.00000000001
    error('The original graph has more than %d connected component', c);
end;

% 假设F固定，
% updateA
% 逐行update A
for iter = 1:NITER
    %F之间的L2距离
    distf = L2_distance_1(F',F');
    A = zeros(num);
    for i=1:num
        if islocal == 1
            %islocal idxa0表示离xi最近的点的坐标
            idxa0 = idx(i,2:k+1);
        else
            idxa0 = 1:num;
        end;
        dfi = distf(i,idxa0);   % 选出与fi相距最近的k个distf距离
        dxi = distX(i,idxa0);   % 选出与xi相距最近的k个distX距离
        ad = -(dxi+lambda*dfi)/(2*r);       % r is fixed ;ad = d_i/(2 r); d_ij = d_ij^x + r * d_ij^f
        A(i,idxa0) = EProjSimplex_new(ad);  % s_ij = -ad + n; --> PCAN:p980 (28)
    end;

    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;
    F_old = F;
    [F, temp, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.00000000001
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
    else
        break;
    end;

end;

%[labv, tem, y] = unique(round(0.1*round(1000*F)),'rows');
[clusternum, y]=graphconncomp(sparse(A)); y = y';
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', c)
end;

