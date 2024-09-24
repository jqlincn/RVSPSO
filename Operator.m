function [V,CanSol] = Operator(Problem,Arc,V0,V,k,beta)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
 
% This function is written by Jianqing Lin
    
    %% Environmental Selection
    NP = Problem.N;
    D  = Problem.D;
    M  = Problem.M;
    theta = Problem.FE/Problem.maxFE;

    [index,Cindex] = KrigingSelectTEST(Arc.objs,V,k,theta);
    RefPop         = Arc(index);

    CanSolDec = [];

    for v = 1 : numel(index)

        RefSol = RefPop(v);
        Vector = RefSol.objs;
        
        CluSol = Arc(Cindex==v);
        [~,FrontNo,~] = EnvironmentalSelection(CluSol,size(CluSol,2));

        GoodSol = CluSol(FrontNo==1);
        BadSol  = CluSol(FrontNo == max(FrontNo));

        %% Statistics Mean Value
        Good_Model_Mean = mean(GoodSol.decs,1);
        Bad_Model_Mean  = mean(BadSol.decs,1);
        Model_Dif       = Good_Model_Mean - Bad_Model_Mean;
        Model_Dif_sort  = sort(abs(Model_Dif),'descend');
        Index_dif  = find(abs(Model_Dif) > Model_Dif_sort(ceil(beta*D)));

        %% Surrogate Optimization
        ArcDec = Arc.decs;
        
        % Variables Selection
        Decs_Surrogate = ArcDec(:,Index_dif);
        Objs_Surrogate = Arc.objs;
        
        % Build RBF-Surrogate 
        para = RBFBuild(Problem,Decs_Surrogate,Objs_Surrogate);
    
        % Reproduction
        [Population,~,~] = EnvironmentalSelection(Arc,NP);
        PopDec           = Population.decs;
        Pop_Surrogate    = Surrogate_individual(PopDec(:,Index_dif),Population.objs);
        Off_Surrogate    = Reproduction(Problem,Pop_Surrogate,RefSol,para,Index_dif);
        
        OffObj = Off_Surrogate.objs;
        OffObj = OffObj - repmat(min(OffObj,[],1),size(OffObj,1),1);
        Angle  = acos(1-pdist2(OffObj,Vector,'cosine'));
        APD    = (1+M*theta*Angle).*sqrt(sum(OffObj.^2,2));

        % Select the one with the minimum APD value
        [~,best] = min(APD);
        OffSol_d = Off_Surrogate(best);

        RefSolDec            = RefSol.decs;
        RefSolDec(Index_dif) = OffSol_d.decs;

        CanSolDec = [CanSolDec;RefSolDec];
        
    end
    
    CanSol  = Problem.Evaluation(CanSolDec);

    % Update V
    V(1:size(V0,1),:) = V0.*repmat(max(CanSol.objs,[],1)-min(CanSol.objs,[],1),size(V0,1),1);
end