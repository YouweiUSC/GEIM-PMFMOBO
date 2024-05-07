function RVEA(Global)
% <algorithm> <R>
% Reference vector guided evolutionary algorithm
% alpha ---   2 --- The parameter controlling the rate of change of penalty
% fr    --- 0.1 --- The frequency of employing reference vector adaptation
% p_driven --- 0 --- 1 for preference-driven and 0 for general

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, M. Olhofer, and B. Sendhoff, A reference vector guided
% evolutionary algorithm for many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2016, 20(5): 773-791.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [alpha,fr,p_driven] = Global.ParameterSet(2,0.1,0);

    %% Generate the reference points and random population
    preference_driven = p_driven;
    if preference_driven
        switch Global.M
            case 2   
                H_ref1 = [3;3;3];
                H_ref2 = 0;
                s2o2=sqrt(2)/2;
                Vc=[0,1;1,0;s2o2,s2o2;];
                Rc=[0.1;0.1;0.1;];   
                n_par = length(Rc);V0=[];
                for i=1:n_par
                    V0 = [V0;GenPerVec(H_ref1(i), H_ref2, Global.M, Vc(i,:), Rc(i))];
                end
                Global.N = length(V0(:,1));
            case 3
                H_ref1 = 1*ones(4,1);
                H_ref2 = 0;
                s2o2=sqrt(3)/3;
                Vc = [0,0,1;1,0,0;0,1,0;s2o2,s2o2,s2o2;];  
                Rc = 0.09*ones(4,1);  
                n_par = length(Rc);V0=[];
                for i=1:n_par
                    V0 = [V0;GenPerVec(H_ref1(i), H_ref2, Global.M, Vc(i,:), Rc(i))];
                end
                Global.N = length(V0(:,1));
        end
    else
        [V0,Global.N] = UniformPoint(Global.N,Global.M);
    end
    V             = V0;
    Population    = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = randi(length(Population),1,Global.N);
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^alpha);
        if ~mod(Global.gen,ceil(fr*Global.maxgen))
            V(1:Global.N,:) = ReferenceVectorAdaptation(Population.objs,V0);
        end
    end
end