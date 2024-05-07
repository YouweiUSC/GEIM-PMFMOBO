H_ref1 = [8;8;8];
H_ref2 = 0;s2o2=sqrt(2)/2;
Vc=[0,1;1,0;s2o2,s2o2;];
Rc=[0.12;0.15;0.12;];   
n_par = length(Rc);V=[];
for i=1:n_par
    V = [V;GenPerVec(H_ref1(i), H_ref2, num_obj, Vc(i,:), Rc(i))];
end