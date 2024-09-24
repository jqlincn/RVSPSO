classdef RVSPSO < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Reference Vector Guided Variables Selection for Expensive Large-Scale Multiobjective Optimization
% k  ---  5  --- The number of re-evaluated solutions

%------------------------------- Reference --------------------------------
% J. Lin, C. He, X. Liu, and L. Pan. Reference Vector Guided Variables 
% Selection for Expensive Large-Scale Multiobjective Optimization. In: 2024 
% IEEE Congress on Evolutionary Computation (CEC), Yokohama, Japan, 2024, 
% pp. 1-8, doi: 10.1109/CEC60901.2024.10611889.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jianqing Lin

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [k] = Algorithm.ParameterSet(5);
            NI       = 100;     % Initial number of solutions
            
            [V0,~] = UniformPoint(Problem.N,Problem.M);
	        V      = V0;

            %% Generate initial population
            InitPopDec = repmat((Problem.upper - Problem.lower),NI,1) .* lhsdesign(NI, Problem.D) + repmat(Problem.lower, NI, 1);
            Arc        = Problem.Evaluation(InitPopDec);
            
            %% Optimization
            while Algorithm.NotTerminated(Arc)
                beta        = 1-((Problem.FE - NI)/(Problem.maxFE - NI))*0.3;
                [V,CanSol]  = Operator(Problem,Arc,V0,V,k,beta);
                Arc         = [Arc,CanSol];
            end
        end
    end
end