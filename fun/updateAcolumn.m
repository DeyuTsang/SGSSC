function [anew] = updateAcolumn(d)
% d: one column in D matrix corresponding to a column in matrix A. d must be all nonnegative
% This function solves the following proble: 
%  min. 1/2|| x - d ||_2^2 st x>=0, sum(x) = 1

if sum(d<0)>1
    error('D must be a nonnegative matrix!');
    return
end
n = length(d);
[d,idx] = sort(d);
k = 0;
for i=1:n
    lambda = mean(d(1:i)) + 1/i;
    if sum(lambda > d(1:i)) == i && sum(d(i+1:end)>lambda)==n-i
        k = i;
        break;
    end;
end;
if k>0
    anew = d*0;
    anew(idx(1:k)) = -d(1:k) + lambda;
else
    error('Something is wrong!')
    return
end;