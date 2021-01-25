function BiCo(Global)
% <algorithm> <H-N>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
   % [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    ArcPop = [];
    
    
    %% Optimization
    while Global.NotTermination(Population)
        
     AllPop = [Population,ArcPop];  
     MatingPool = MatingSelection2(Population,ArcPop,Global.N);
     Offspring  = Global.Variation(MatingPool(1:Global.N)); 
     ArcPop = UpdateArc([AllPop,Offspring],Global.N); 
     Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
    
end