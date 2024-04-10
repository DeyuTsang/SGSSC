function [x ft] = EProj(v)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. 1>= x >=0

% Ming Yin
% 2016-04

n = length(v);

for i =1:n
   if v(i)>1
       x(i) = 1;
   elseif v(i)<0
        x(i) = 0;
   else
        x(i) = v(i);
   end       
end

