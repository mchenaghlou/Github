function [ indm, centers, NIs, uq1, Radiis] = KMeansStep(curr_ob, centers, NIs, Radiis)
%KMEANSSTEP Summary of this function goes here
%   Detailed explanation goes here

uq1 = calculate_cluster_membership(curr_ob, centers, 'binary');
alld=pdist2(curr_ob,centers);
[m,indm]=min(alld);
NIs(indm)=NIs(indm)+1;
Radiis(indm) = Radiis(indm) + (1/NIs(indm))*m;
centers(indm,:)=centers(indm,:)+(1/NIs(indm))*(curr_ob-centers(indm,:));

end

