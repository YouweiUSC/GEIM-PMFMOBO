classdef ComprDes < PROBLEM
% <multi> <real> <large/none> <expensive/none>
% Compressor design problem based on Multall and FFD

%------------------------------- Reference --------------------------------
% E. Zitzler, K. Deb, and L. Thiele, Comparison of multiobjective
% evolutionary algorithms: Empirical results, Evolutionary computation,
% 2000, 8(2): 173-195.
%------------------------------- Copyright --------------------------------
    properties
        GeoPar
        Profile
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 2;
            if isempty(obj.D); obj.D = 144; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = 50*ones(1,obj.D);
            obj.encoding = 'real';
            % define the FFD box
            [obj.GeoPar,obj.Profile] = FFDBoxDef(6);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec,varargin)
            nPop = size(PopDec,1);fidelity = 6*ones(nPop,1);
            AllData=RunMultall(PopDec,fidelity,obj.GeoPar,obj.Profile);
            PopObj = -1.*[AllData.pressure_ratio' AllData.poly_eff_tt'];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R(:,1) = linspace(0,1,10)';
            R(:,2) = 1 - R(:,1).^0.5;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = obj.GetOptimum(100);
        end
    end
end