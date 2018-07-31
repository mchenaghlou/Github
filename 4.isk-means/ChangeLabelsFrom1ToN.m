function [ labels ] = ChangeLabelsFrom1ToN( labels )
%CHANGELABELSFROM1TON Summary of this function goes here
%   Detailed explanation goes here

u_labels = unique(labels, 'stable');
i = 1;
while i <= length(u_labels)
    if i == u_labels(i)
        i = i + 1;
        continue;
    else
        labels(find(labels == i)) = 0;
        labels(find(labels == u_labels(i))) = i;
        labels(find(labels == 0)) = u_labels(i);
    end
    u_labels = unique(labels, 'stable');
    i = i + 1;
end

end

