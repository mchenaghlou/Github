%index
function [ XB1,Mindist,GVec,CVec,PsiVec] = IXB(x,u1,oCenters,Centers,lambda,...
    GVec,CVec,PsiVec,k, N,oMindist)
% x: datapoint
% u1:
% oCenters:
% Centers:
% lambda:
% GVec:
% CVec:
% PsiVec:
% k:
% N:
% oMindist:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PE(fixindex+i) = sum(sum(U.*log(U)))/(-1*n*log(k));
% sum(sum(U.^2.*pdist2(centers,data).^2))/n/min(min(pdist(centers).^2));
% N=N+1;
% u1=zeros(k,1);
% u1(ind)=1;

% this is h_{q}
Mindist = min(min(pdist(Centers).^2));
if(numel(Mindist)==0 || Mindist==0)
    [~,ind] = max(u1);
    Mindist= max([oMindist norm(x-Centers(ind,:)).^2]);
end


try
    for i=1:k
        [CVec(i),GVec(i,:),PsiVec(i)] = iCompact(x,u1(i),oCenters(i,:),Centers(i,:),lambda,CVec(i),GVec(i,:), PsiVec(i));
    end
catch E
    milad = 1;
end
if(lambda<1)
    XB1 = ((1-lambda)*sum(CVec))/Mindist;
else
    XB1 = sum(CVec)/(N*Mindist);
end

