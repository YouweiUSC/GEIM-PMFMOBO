function U=UpdateLUBs(U,p)
% INPUT
% U: the LUBs to be updated
% p: the new Pareto point
% OUTPUT
% U: the updated LUBs
% -----------------------------
% n: number of LUBs, M: number of sub-objectives
[m,M]=size(U);
D = max(repmat(p,m,1)-U,[],2);
% current LUBs strongly dominated by p
SDU = U(D<0,:);
% current LUBs dominated by p but nor strongly
WDU = U(D==0,:);
% delete the LUBs strongly dominated by p
U(D<0,:) = [];
% take the projection points of p on the boundary of
% each  [-Inf, SDU(i,:)] as the candidates of new LUBs
[n,~] = size(SDU);
C = repmat(SDU,M,1);
for i=1:M
    C(n*(i-1)+1:n*i,i)=repmat(p(i),n,1);
end
% Update the LUBs
% U = [U;-ParetoSet(-[C;WDU])];% original code
PopObj = -[C;WDU];
index = NDSort(PopObj,1)==1;
U = [U;-PopObj(index',:)];% PlatEMO-based implementation
end