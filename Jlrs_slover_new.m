function [R,E,B, err, obj] =  Jlrs_slover(X, k, Z1, Z2, Z3, beta, rho, Iter, method, Debug )
% Jointly learning representation and similarity of data, optimized by LADMSAP
%
%Input: 
%    X :  data with size of m x n; m is the dimension, while n is the
%    number of samples. 
%    k :   the number of the nearest neighbors
%    Z1, Z2, Z3: regularized para. 
%    Iter: the number of interations
%    method: ssc or lrr;
%   Debug: display the err and change for each iteration
%
%Output: 
%   R: representation matrix;
%   B: affinity matrix;
%   E: the residual error
%
% Ming Yin
% 2016-02


if (~exist('k','var'))
    k = 10;
end

if (~exist('Iter','var'))
    Iter = 200;
end

if (~exist('method','var'))
    method = 'ssc'; % or 'lrr'
end

if (~exist('Debug','var'))
    Debug = 1;
end

if (~exist('Z1','var'))
   Z1 = 100;
end

if (~exist('Z2','var'))
   Z2 = 0.01;
end

if (~exist('Z3','var'))
    Z3 = 0.001;
end

if (~exist('rho','var'))
    rho = 1.1;
end

[m,n] = size(X);

%-------------------  init  ---------------------
ro  = rho;
beta_max = 1e5;
if (~exist('beta','var')|| isempty(beta))
    beta = norm(X)*min(m,n)*1e-4/10;
end
% beta = norm(X)*min(m,n)*1e-4/10;

lambda = zeros(n+m,n);
edta_E =  1;
edta_R =  4*1.01*norm(X).^2;      %1.02 * 3 * (norm(X)^2);
edta_Q =  1;       
tol1= 1e-4;
tol2 = 1e-3;
 
A_R = [X;eye(n,n)];
A_E = [eye(m,m);zeros(n,m)];
A_Q = [zeros(m,n);-eye(n,n)];
b   = [X;zeros(n,n)];

err  = [];
obj  = [];
E = zeros(m,n);
R = zeros(n,n);
Q = R;
B = R;
LB = diag(sum(B,2)) - (abs(B)+abs(B'))/2;
idx_adjust = repmat( [0:n-1].*n,k+1,1);
norm_b = norm(b,'fro');
 for t =  1:Iter

        Et = E;  Rt = R;  Qt = Q;   Bt = B;
        %------------ step 1
        sumAx = A_R*R + A_E*E + A_Q*Q ;
        lambda_hat = lambda + beta*(sumAx - b);        

        %----------- step 2
        %--------- E
        sigma_E = edta_E * beta; 
        W_E = E - (A_E' * lambda_hat) ./ sigma_E;
        E = ( sigma_E / (2+sigma_E) ) .* W_E;

        %--------- R        
        sigma_R = edta_R*beta;
        W_R       = R - (A_R'*lambda_hat)./sigma_R;   %(edta_R*beta);
        
        if strcmp(method,  'ssc')  
            R = sign(W_R) .* max( abs(W_R) - Z1/sigma_R,0);
            R = R - diag(diag(R));
        else
        % nuclear-norm regularzier,  'lrr'
            p = Z1/sigma_R;        
            R =  subfun_update_Mlowrank(W_R, 1/p );
            
%             [U,S,V] = svd(W_R,'econ');
%             S = diag(S);
%             svp = length(find(S> p));
%             if svp>=1
%                 S = S(1:svp) - p;
%             else
%                 svp = 1;
%                 S = 0;
%             end
%             R = U(:,1:svp)*diag(S)*V(:,1:svp)';   
            
        end

        %---------- Q
        sigma_Q = edta_Q*beta / Z2;
        W_Q = Q - A_Q' * lambda_hat / (edta_Q*beta);
        Q = (LB+LB' + (sigma_Q * eye(size(LB)) ) );
        Q = Q\(sigma_Q * W_Q) ;

        %---------- B
%         B = zeros(n,n);
%         for i = 1:n 
%             di = (Z2 /(4 * Z3))*sum((repmat(Q(:,i),1,n) - Q).^2,1)';
%             [sorted_di,idx] = sort(di);
%             B(idx(2:k+1),i) = ( (1 + sum( sorted_di(2:k+1) ) ) / k )* ones(k,1) - di(idx(2:k+1));
%         end
%         B = max( B,0 );
%         B = (B + B') / 2;
        
        B = zeros(n,n);
        [knnidx, knndist] =  knnsearch(Q',Q','dist','euclidean','K',k+1);
        knnidx = knnidx';knndist =knndist';
        knnidx = knnidx + idx_adjust;
        dQ =  (Z2/4/Z3).*(knndist.^2);
        B (knnidx(2:k+1,:)) = repmat((1 + sum(dQ(2:k+1,:),1))./k,k,1) - dQ(2:k+1,:);
        B = max( B  ,0 );
        % balance B 
        B  = (B  + B') / 2;

        sumAx = A_R*R + A_E*E + A_Q*Q ;        
        err1     = norm(sumAx - b,'fro')/norm_b;
        err       = [err err1 ];
        
%      regch = max([ sqrt(edta_E) * norm(E - Et,'fro'), sqrt(edta_Q) * norm(Q - Qt,'fro'),...
%      sqrt(edta_R) * norm(R - Rt,'fro') ])/norm(b,'fro');     %beta * , norm(B - Bt,'fro')    
       regch = max([  norm(E - Et,'fro'),  norm(Q - Qt,'fro'), norm(R - Rt,'fro'), norm(B - Bt,'fro')  ])/norm_b;     %beta * 
              
        if Debug
           fprintf('ladmpsap: iter:%d,    err1:%f,	 regch :%f,   beta = %f \n\r', t,  err1,  regch,  beta);
           %rank(R)= %d,  ro = %f \n\r', t,  err1,  regch,  rank(R),  ro);
        end
                
        LB     = diag(sum(B,2)) - (abs(B)+abs(B'))/2;
        obj(t) = Z1*(norm(R,'fro').^2+ trace(R*LB*R')+norm(B,'fro').^2 )+ norm(E,'fro').^2;
       
       if err1< tol1 && beta*regch <tol2
            break;
       end
       
       %------- updating lambda beta
       lambda = lambda + beta *(sumAx - b);
       if (beta*regch < tol2)    
            ro = rho;
       else
           ro = 1;
       end
       beta      = min(beta_max, ro*beta);   
 end
 
 function [J ] = subfun_update_Mlowrank(X, mu )
% this functoin solves the  following prob:
% J = argmin
%        || J ||_* + 0.5*mu*||J - ( X )||_F^2

% if ~silent
%     fprintf('sin val >  thre =%.3f NORMAL SVD', 1/mu );
% end

[U,sigma,V] = svd(   X,'econ');

sigma  = diag(sigma)';
ind       = sigma>1/mu;

if  sum(ind)~=0
    sigma = sigma( ind )-1/mu;
    J         = bsxfun(@times, U(:,ind) , sigma(ind) ) * ( V(:,ind)') ;
else
    J = zeros( size( X ),'single');
end

 
   