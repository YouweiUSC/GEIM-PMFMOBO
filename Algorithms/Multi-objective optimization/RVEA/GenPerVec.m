function V = GenPerVec(H_ref1, H_ref2, num_obj, Vc, Rc)
[N,Vs] = F_weight(H_ref1, H_ref2, num_obj); % N ~ number of reference vectors
for i = 1:N
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end
V = Rc*Vs+(1-Rc)*Vc;
for i = 1:N
    V(i,:) = V(i,:)./norm(V(i,:));
end
end

%Function to generate uniformly distributed weight vectors
function [N,W] = F_weight(p1,p2,M)
    [N,W] = T_weight(p1,M);
    if p2 > 0
        [N2,W2] = T_weight(p2,M);
        N = N+N2;
        W = [W;W2*0.5+(1 - 0.5)/(M)];
    end
end

function [N,W] = T_weight(H,M)
    N = nchoosek(H+M-1,M-1);
    Temp = nchoosek(1:H+M-1,M-1)-repmat(0:M-2,nchoosek(H+M-1,M-1),1)-1;
    W = zeros(N,M);
    W(:,1) = Temp(:,1)-0;
    for i = 2 : M-1
        W(:,i) = Temp(:,i)-Temp(:,i-1);
    end
    W(:,end) = H-Temp(:,end);
    W = W/H;
end

