function ArcPop = UpdateArc(Population,N)
 

 [FrontNo,MaxFNo] = NDSort([Population.objs,sum(max(0,Population.cons),2)],1);
  Temp1 = FrontNo==1;
  
%   [~,index]=sort(sum(max(0,Population.cons),2));
%   Temp2 = find(index<N+1);
%   Temp2 = index(Temp2);
%   
%   Temp1(Temp2')=false;
  
  Population = Population(Temp1==1);
  
  Temp = find(sum(max(0,Population.cons),2)>0);
  Population = Population(Temp);
  
  if length(Population)<N 
      ArcPop = Population;
  else
      
    Zmax     = max(Population.objs,[],1);
    Next(1:size(Population,2)) = true;
    % Select the solutions in the last front
    Delete = LastSelection(Population(Next).objs,-sum(max(0,Population.cons),2),sum(Next)-N,Zmax);
    Temp = find(Next);
    Next(Temp(Delete)) = false;
     ArcPop = Population(Next);
      
     % ArcPop = Angleslection(Population,N);
  end
  
%   ArcFrontNo = [];
%   ArcCrowdDis = [];
%   if length(ArcPop)~= 0
%       %ArcFrontNo = [1:length(ArcPop)];
%       ArcFrontNo = ones(1,length(ArcPop));
%       ArcCrowdDis =ones(1,length(ArcPop));% [1:length(ArcPop)];
%       %ArcFrontNo = (N+1)*ones(1,length(ArcPop));  
%      % ArcCrowdDis = CrowdingDistance(ArcPop.objs,ArcFrontNo);
%      % ArcCrowdDis =  [N+1:N+length(ArcPop)];
%   end  
  %Next1 = find(Population(Next==1).cons>0);
end


function Delete = LastSelection(PopObj,PopCons,K,Zmax)
% Select part of the solutions in the last front
    [N,M]  = size(PopObj);
    
    PopObj = (PopObj-repmat(Zmax,N,1))./(repmat(min(PopObj),N,1)-repmat(Zmax,N,1)- 1e-10);
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
   % PopObj = real(PopObj);
      Cosine   = 1 - pdist2(PopObj,PopObj,'cosine');
      Cosine   = Cosine.*(1-eye(size(PopObj,1)));
%       SDE = zeros(1,N);
%       for i=1:N
%         SPopuObj = PopObj;
%         Temp     = repmat(PopObj(i,:),N,1);
%         Shifted  = PopObj < Temp;
%         SPopuObj(Shifted) = Temp(Shifted);
%         Distance = pdist2(PopObj(i,:),SPopuObj);
%         [~,index] = sort(Distance,2);
%         SDE(i) = Distance(index(floor(sqrt(N))+1));
%       end
          
      
    
    %% Environmental selection
    Delete  = false(1,N);
    % Select K solutions one by one
    while sum(Delete) < K
        [Jmin_row,Jmin_column] = find(Cosine==max(max(Cosine)));
        j = randi(length(Jmin_row));
        Temp_1 = Jmin_row(j);
        Temp_2 = Jmin_column(j);
         
        if  (PopCons(Temp_1)<PopCons(Temp_2)) ||(PopCons(Temp_1)==PopCons(Temp_2) && rand<0.5)
            Delete(Temp_1) = true;
            Cosine(:,Temp_1)=0;
            Cosine(Temp_1,:)=0;
        else
            Delete(Temp_2) = true;
            Cosine(:,Temp_2)=0;
            Cosine(Temp_2,:)=0;
        end
    end
end
% 
% function ArcPop = Angleslection(Population,N)
% 
% end