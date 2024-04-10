
% transform a label information into a binary weight matrix
% Ming Yin,
 
function A = lable2weight(trgnd,n)

A = eye(n);
c = unique(trgnd);  
for i = 1:length(c)
    idx = find(trgnd==c(i));
    A(idx,idx) = 1;
end