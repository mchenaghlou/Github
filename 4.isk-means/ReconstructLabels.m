function [ c,s ] = ReconstructLabels( Clusters,  S)
%RECONSTRUCTLABELS Summary of this function goes here
%   Detailed explanation goes here

X = Clusters;
hi = histcounts(X, numel(unique(X)));
singles = find(hi==1);
inds = ismember(X,singles);
X(find(inds)) = 0;

ss = unique(X);
i = 0;
while i < length(ss)-1
    i = i+1;
    if ismember(i,X)
        continue
    else
        j = ss(i+1);
        t = find(X==j);
        X(t) = i;
        ss = unique(X);
        i = i -1;
    end
end
c = X';
t = unique(c);
c = c';

X =S;
ss = unique(X);
i = 1;
if ismember(0,X)
    count = length(ss) -1;
else
   count = length(ss);
end
while i <= count
    if ismember(i,X)
        i = i + 1;
        continue
    else
        if i == length(ss)
            j = ss(i);
        else
            j = ss(i+1);
        end
        
        t = find(X==j);
        X(t) = i;
        ss = unique(X);
        i = i -1;
    end
    i = i+1;
end
s = X';
t = unique(s);
s= s';


% c = Clusters;
% t=unique(c);
% 
% for i = 1:max(t)
%     if ismember(i, t) == false
%         c = sortClusterLabels(c, i, t);
%     end
% end
% 
% tc=unique(c);
% k=numel(tc);
% c(c==0) = k;
% % t=unique(c);
% 
% s = S;
% ts=unique(s);
% for i = 1:max(ts)
%     if ismember(i, ts) == false
%         s = sortClusterLabels(s, i, ts);
%     end
% end
% 
% ts=unique(s);
% k=numel(ts);
% s(s==0) = k;


end

