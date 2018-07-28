function [XBLambdaq1, hq1, uq1] = IXB_ff_milad( xq1, Vq1, hq, XBLambdaq, lambda, k, type)
%IXB_MILAD Summary of this function goes here
%   Detailed explanation goes here

hq1 = calculate_Hq(Vq1, xq1, hq);
uq1 = calculate_cluster_membership(xq1, Vq1, type);

DeltaJq1 = 0;
for i = 1:k
    DeltaJq1 = DeltaJq1  + (uq1(i)^2) * (pdist2(xq1, Vq1(i, :))^2);
end
XBLambdaq1 = (1/hq1) * (lambda * hq * XBLambdaq + (1-lambda)*DeltaJq1);
% XBLambdaq1 =   (1/1)   * (lambda * 1 * XBLambdaq + (1-lambda)*DeltaJq1);
end