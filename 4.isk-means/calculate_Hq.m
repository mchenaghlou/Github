function [ hq ] = calculate_Hq(V, xq1, hq_prev)
%CALCULATE_HQ Summary of this function goes here
%   Detailed explanation goes here

if(size(V, 1) == 1)
   [max_val, max_ind] = max(pdist2(xq1, V)^2); 
   hq = max(max_val, hq_prev);
else
   hq = min(min(pdist(V).^2));  
end


end


