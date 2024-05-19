classdef PGEIM_MOBO < ALGORITHM
% <multi> <real> <expensive>
% PGEIM_MOBO-EGO: parallel generalized expected improvement matrix multi-objective bayesian optimization 
% q    ---   6 --- Batch size  
% NewOpt --- 1 --- Start new optimization
% maxIter --- 100 --- max interation for the inner optimizer
% HighDim --- 0 --- High-dimension problem indicator
% CR --- 6 --- Cost ratio

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [q,NewOpt,maxIter,HighDim,CR] = Algorithm.ParameterSet(6,1,100,0,6);
            switch class(Problem)
                case 'ComprDes'
                    NewOpt = 0;
                otherwise
                    NewOpt = NewOpt;
            end
            % number of initial design points
            switch Problem.D
                case 6
                    n_LF = 72;
                    n_HF = 20;
                otherwise 
                    n_LF = 10*Problem.D;
                    n_HF = 4*Problem.D;
            end
           %% Generate the initial design points
            if NewOpt
                PopDecHF = Problem.lower + (Problem.upper-Problem.lower).*lhsdesign(n_HF,Problem.D,'criterion','maximin','iterations',1000);
                PopDecLF = Problem.lower + (Problem.upper-Problem.lower).*lhsdesign(n_LF,Problem.D,'criterion','maximin','iterations',1000);
                PopulationHF = Problem.Evaluation(PopDecHF,1);  
                PopulationLF = Problem.Evaluation(PopDecLF,CR);  
                [FrontNo,~]  = NDSort(PopulationHF.objs,1); % find non-dominated solutions (FrontNo==1)
                Algorithm.NotTerminated(PopulationHF(FrontNo==1));  % store pupulation 
                % 30 for ZDT1-3, 28 for DTLZ5
                switch class(Problem)
                    case {'ZDT1','ZDT2','ZDT3'}
                        N_iter = ceil(9*30/q);
                    case 'DTLZ5'
                        N_iter = ceil(9*28/q);
                    otherwise
                        N_iter = 30;
                end
            else
                % ########################### for compressor optimization
                load('test_x.mat','test_x')
                load('test_y.mat','PR','Eff')
                Ind = randperm(n_LF,n_HF)';
                PopDecHF = test_x(Ind,:);
                PopDecLF = test_x;
                PopObjHF = -1.*[PR(6,Ind)' Eff(6,Ind)'];
                PopObjLF = -1.*[PR(3,:)' Eff(3,:)'];
                clear test_x PR Eff
                N_iter = 30;
            end
            NCurrentSample = n_LF+n_HF;
            AllowFE  = N_iter*q+NCurrentSample;
            StopCri  = NCurrentSample<AllowFE;
            %% Optimization
            IterCount = 0;
            while StopCri
                IterCount = IterCount+1;
                fprintf('No. Iteration: %d\n',IterCount)
              %% Scale the objective values 
                if NewOpt
                    [PopDecHF,PopObjHF,PopDecLF,PopObjLF] = ...
                        deal(PopulationHF.decs,PopulationHF.objs,PopulationLF.decs,PopulationLF.objs); 
                end
                PopObjHFScaled = PopObjHF; 
                PopObjLFScaled = PopObjLF;                 
              %% Bulid HK model for each objective function 
                % bulid Hierarchical Kriging models for all the objective functions
                disp('     build model')
                sample_x = {PopDecLF;PopDecHF};
                HKmodel  = cell(1,Problem.M);
                optionsHK.hyperest='Boxmin'; % optionsHK.hyperest='';  Boxmin DIC
                if HighDim
                    optionsHK.hyperest='DIC';
                end
                for i = 1 : Problem.M
                    sample_y = {PopObjLFScaled(:,i);PopObjHFScaled(:,i)};
                    HKmodel{i}= train_MLHK(sample_x,sample_y, optionsHK);
                end
            %% determine the idea and nadir points
                [Pidea, Pnadir]= IdealNadir(Problem,HKmodel);
                PopObjHFScaled = (PopObjHF-Pidea)./(Pnadir-Pidea); 
                PopObjLFScaled = (PopObjLF-Pidea)./(Pnadir-Pidea); 
                index          = NDSort(PopObjHFScaled,1)==1;
                NDS_scaled     = PopObjHFScaled(index',:);
                for i = 1 : Problem.M
                    sample_y = {PopObjLFScaled(:,i);PopObjHFScaled(:,i)};
                    HKmodel{i}= train_MLHK(sample_x,sample_y, optionsHK);
                end
            %%  Maximize PGEIM and select q candidate points
              disp('     Obtain infill samples')   
              SelectDecs = Opt_PGEIM(Problem,HKmodel,NDS_scaled,maxIter,q); 
            %% Aggregate data
              disp('     Evaluate samples')
                IndividualNew = Problem.Evaluation(SelectDecs);
                if NewOpt
                    PopulationHF = [PopulationHF,IndividualNew];
                    [FrontNo,~] = NDSort(PopulationHF.objs,1); % find non-dominated solutions 
                    Algorithm.NotTerminated(PopulationHF(FrontNo==1)); 
                else
                    PopDecHF = [PopDecHF;SelectDecs];
                    PopObjHF = [PopObjHF;IndividualNew.objs];                 
                    save(['Data\Res_maxFE_' mfilename '_' num2str(Problem.maxFE) '_' datestr(datetime('now'),'mm-dd-HH-MM')], ...
                        'PopDecHF','PopObjHF','PopDecLF','PopObjLF');
                end
                NCurrentSample= NCurrentSample+q;
                StopCri = NCurrentSample<AllowFE;
            end
            
        end
    end
end
