classdef OSY < PROBLEM
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
            obj.lower    = [0 0 1 0 1 0];
            obj.upper    = [10 10 5 6 5 10];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec,varargin)
            x1 = PopDec(:,1);
            x2 = PopDec(:,2);
            x3 = PopDec(:,3);
            x4 = PopDec(:,4);
            x5 = PopDec(:,5);
            x6 = PopDec(:,6);
            if isempty(varargin)
                PopObj(:,1) = -25*(x1-2).^2-(x2-2).^2-(x3-1).^2-(x4-4).^2-(x5-1).^2;
                PopObj(:,2) = x1.^2+x2.^2+x3.^2+x4.^2+x5.^2+x6.^2;
            else
                % for multi-fidelity cases
                CostRatio=varargin{1};
                if CostRatio>1% LF
                    PopObj(:,1) = -20*(x1-2).^2-2*(x2-1).^2-3*(x3-1).^2-(x4-4).^2-(x5-1).^2;
                    PopObj(:,2) = 0.9*x1.^2+0.8*x2.^2+1.1*x3.^2+1.2*x4.^2+x5.^2+x6.^2+x1.*x2-x5;
                else% HF
                    PopObj(:,1) = -25*(x1-2).^2-(x2-2).^2-(x3-1).^2-(x4-4).^2-(x5-1).^2;
                    PopObj(:,2) = x1.^2+x2.^2+x3.^2+x4.^2+x5.^2+x6.^2;
                end
            end
        end
         %% Calculate Constraint values
        function PopCon = CalCon(obj,PopDec,varargin)
            x1 = PopDec(:,1);
            x2 = PopDec(:,2);
            x3 = PopDec(:,3);
            x4 = PopDec(:,4);
            x5 = PopDec(:,5);
            x6 = PopDec(:,6);
            if isempty(varargin)
                PopCon(:,1) = (x3-3).^2+(x4-4);
                PopCon(:,2) = -(x5-3).^2-x6+4;
                PopCon(:,3) = x1-3*x2-2;
            else
                % for multi-fidelity cases
                CostRatio=varargin{1};
                if CostRatio>1% LF
                    PopCon(:,1) = (x3-2).^2+2*(x4);
                    PopCon(:,2) = -(x5-1).^2-x6;
                    PopCon(:,3) = x1-2*x2-1;
                else% HF
                    PopCon(:,1) = (x3-3).^2+(x4-4);
                    PopCon(:,2) = -(x5-3).^2-x6+4;
                    PopCon(:,3) = x1-3*x2-2;
                end
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [0 280];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end