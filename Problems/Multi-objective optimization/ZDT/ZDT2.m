classdef ZDT2 < PROBLEM
% <multi> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Zitzler, Deb, and Thiele

%------------------------------- Reference --------------------------------
% E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 6; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec,varargin)
            if isempty(varargin)
                PopObj(:,1) = PopDec(:,1);
                g = 1 + 9*mean(PopDec(:,2:end),2);
                h = 1 - (PopObj(:,1)./g).^2;
                PopObj(:,2) = g.*h;
            else
                % for multi-fidelity cases
                CostRatio=varargin{1};
                if CostRatio>1% LF
                    PopObj(:,1) = 0.8*PopDec(:,1)+0.5;
                    g = 1 + 9*mean(PopDec(:,2:end),2);
                    h = 1 - (PopObj(:,1)./g).^2;
                    PopObj(:,2) = (0.9*g+0.2).*(1.1*h-0.2)+2*PopDec(:,1);
                else% HF
                    PopObj(:,1) = PopDec(:,1);
                    g = 1 + 9*mean(PopDec(:,2:end),2);
                    h = 1 - (PopObj(:,1)./g).^2;
                    PopObj(:,2) = g.*h;
                end
            end
        end
       %% Calculate Constraint values
        function PopCon = CalCon(obj,PopDec,varargin)
            PopCon = zeros(size(PopDec,1),1);
        end
       %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [1.2 1.2];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end