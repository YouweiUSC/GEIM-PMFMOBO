classdef POL < PROBLEM
% <multi> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Zitzler, Deb, and Thiele

%------------------------------- Reference --------------------------------
% E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function Setting(obj)
            obj.M = 2;
            obj.D = 2;
            obj.lower    = -pi*ones(1,obj.Global.D);
            obj.upper    =  pi*ones(1,obj.Global.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec,varargin)
            x = PopDec;
            PopObj(:,2) = (x(:,1)+3).^2 + (x(:,2)+1).^2;
            A1 = 0.5*sin(1)-2*cos(1)+1*sin(2)-1.5*cos(2);
            A2 = 1.5*sin(1)-1*cos(1)+2*sin(2)-0.5*cos(2);
            B1 = 0.5*sin(x(:,1))-2*cos(x(:,1))+1*sin(x(:,2))-1.5*cos(x(:,2));
            B2 = 1.5*sin(x(:,1))-1*cos(x(:,1))+2*sin(x(:,2))-0.5*cos(x(:,2));
            if isempty(varargin)
                PopObj(:,1) = 1+(A1-B1).^2+(A2-B2).^2;
            else
                % for multi-fidelity cases
                CostRatio=varargin{1};
                if CostRatio>1% LF
                    PopObj(:,1) = 1+(0.9*A1-1.2*B1).^2+0.9*(1.2*A2-0.9*B2).^2;
                else% HF
                    PopObj(:,1) = 1+(A1-B1).^2+(A2-B2).^2;
                end
            end
        end
        %% Sample reference points on Pareto front
        function P = GetPF(obj,N)
            load('PF_POL.mat', 'PF_POL')
            P = PF_POL;
        end
    end
end