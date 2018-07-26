function [ initXB, U_nomralized ] = XB_batch(X, V)
%BXB Summary of this function goes here
%   Detailed explanation goes here
n = size(X, 1);
k = size(V, 1);
% initXB = 0;

U_nomralized = zeros(n, k);
for i = 1:size(X, 1)
    U_nomralized(i, :) = calculate_cluster_membership(X(i, :), V, 'fuzzy');
end
dists = pdist2(X, V);
J = 0;
for i = 1:1:n
    for j = 1:1:k
        J = J + (U_nomralized(i, j)^2) * (dists(i, j)^2);
    end
end


if size(V, 1) == 1
    initXB = J / (n * max(dists));
else
    denom = n * min(pdist(V).^2);
    initXB = J / denom;
end


end

