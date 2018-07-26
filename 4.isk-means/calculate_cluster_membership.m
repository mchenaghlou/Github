function [ res ] = calculate_cluster_membership( obs, centers, type)
%CALCULATE_CLUSTER_MEMBERSHIP Summary of this function goes here
%   Detailed explanation goes here
res = zeros(1, size(centers, 1));


if strcmp(type, 'fuzzy') == 1
%     res =  ((alld ./ sum(alld)));
%     if find(res == 0) > 0
%         res2 = zeros(1, length(res));
%         res2(find(res == 0)) = 1;
% 
%     else
%         res = 1 ./ res;
%         res = res / sum(res);
%     end

% This is from Xie-Beni paper.
alld= 1 ./ (pdist2(obs,centers) .^ 2);
if sum(isinf(alld)) == true || sum(isnan(alld) == true) > 0
    res = zeros(1, size(centers, 1));
    [~, index] = ismember(obs,centers,'rows');
    res(index) = 1;
else
    res = alld ./ sum(alld);
end


    % According to the paper: Validity measure for fuzzy clustering by Xie-Beni
%     res = sqrt(1 ./ (alld .^ 2) ) ./ sum(sqrt(1 ./ (alld .^ 2)));
    
    
elseif strcmp(type, 'binary') == 1
    alld= pdist2(obs,centers);
    [~,indm]=min(alld);
    if(length(indm) > 1)
        milad = 1;
    end
    res(indm) = 1;
else
    res = "please select from 'fuzzy' and 'binary'";
end

end

