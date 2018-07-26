% index1
function [C,G,Psi] = iCompact(x,u,Center,NewC,lambda,oC,oG, Psi)
if lambda <1
    T1 = lambda*oC;
    T2 = 2*lambda*(Center-NewC)*oG';
    T3 = lambda*(Center-NewC)*(Center-NewC)'*Psi;
    T4 = u^2*(x-NewC)*(x-NewC)';

    C = T1+T2+T3+T4;

    T1 = lambda*oG;
    T2 = lambda*(Center-NewC)*Psi;
    T3 = u^2*(x-NewC);

    G = T1+T2+T3;
    Psi = lambda*Psi+u^2;
else
    T1 = oC;
    T2 = 2*(Center-NewC)*oG';
    T3 = (Center-NewC)*(Center-NewC)'*Psi;
    T4 = u^2*(x-NewC)*(x-NewC)';

    % T1 = C_{i,n} 
    % T2 = 2*Q_{i,n+1} 
    % T3 = M_{i,n}*B_{i,n+1}  
    % T4 = A_{i, n+1}
    C = T1+T2+T3+T4;

    T1 = oG;
    T2 = (Center-NewC)*Psi;
    T3 = u^2*(x-NewC);

    G = T1+T2+T3;
    Psi = Psi+u^2;
end
end