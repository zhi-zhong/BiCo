function Population = EnvironmentalSelection(Population,N)
% The environmental selection of NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
   % CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    
   % [~,Rank] = sort(CrowdDis(Last),'descend');
   %Next(Last(Rank(1:N-sum(Next)))) = true;
     if sum(Next)+size(Last,2)-N == 0
         Next(Last)=1;
     else
     Del  = Truncation(Population,Last,sum(Next)+size(Last,2)-N);
     Next(Last(~Del)) = true; 
     end
    % Temp = find(Last);

     

    %% Population for next generation
    Population = Population(Next);
   % FrontNo    = FrontNo(Next);
   % CrowdDis   = CrowdDis(Next);
end


function Del = Truncation(Population,Last,K)
% Select part of the solutions by truncation

    

    %% Truncation
    Zmin       = min(Population.objs,[],1);
   % Zmax = max(Population.objs,[],1);
    PopObj = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(Population.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;

%     

    PopObj = (PopObj(Last,:)); 
    
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
      

end