function [affinity] = KNN_CalculateAffinity(data,k)

% set the parameters
sigma = 1;
[n,m] = size(data);
affinity = zeros(n,n);
for i=1:n
    dist = zeros(1,n);
    for j=1:n
       dist(j) = norm(data(i,:)-data(j,:),2)^2;
    end
    [sorted_dist, idx] = sort(dist);
    
    affinity(i,idx(1:k)) = exp(-sorted_dist(1:k)/(2*sigma^2));
   
end




