function U=GenerateLUBs(P,r)
% INPUT
% P: A group of Pareto points
% r: Reference point
% OUTPUT
% U: LUBs
% -----------------------------
% LUB with no Pareto points
U = r;
% Update the LUBs with each Pareto point
[n,~] = size(P);
for i=1:n
    U = UpdateLUBs(U,P(i,:));
end
end