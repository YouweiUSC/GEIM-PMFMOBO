classdef DTLZ5 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 6; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec,varargin)
            g = sum((PopDec(:,2:end)-0.5).^2,2);
            if isempty(varargin)
                PopObj(:,1) = (1+g).*cos(0.5*pi*PopDec(:,1));
                PopObj(:,2) = (1+g).*sin(0.5*pi*PopDec(:,1));
            else
                % for multi-fidelity cases
                CostRatio=varargin{1};
                if CostRatio>1 % LF
                    PopObj(:,1) = (1+0.8*g).*cos(0.5*pi*PopDec(:,1));
                    PopObj(:,2) = (1+1.1*g).*sin(0.5*pi*PopDec(:,1));
                else % HF
                    PopObj(:,1) = (1+g).*cos(0.5*pi*PopDec(:,1));
                    PopObj(:,2) = (1+g).*sin(0.5*pi*PopDec(:,1));
                end
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [1.2 1.2];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M < 4
                R = obj.GetOptimum(100);
            else
                R = [];
            end
        end
    end
end