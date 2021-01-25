function MatingPool = MatingSelection1(Population,ArcPop,N)
% SelectedIndex = zeros(1,N);
  MatingPool = [];
if length(ArcPop)<N
    SelectedIndex = TournamentSelection(2,N,-sum(max(0,Population.cons),2));
    MatingPool = Population(SelectedIndex);
else
    AllPop = [ Population,ArcPop]; 
    Zmin       = min(AllPop.objs,[],1);
   % Zmax = max(Population.objs,[],1);
    PopObj1 = (Population.objs-repmat(Zmin,length(Population.objs),1))./(repmat(max(AllPop.objs),length(Population.objs),1)-repmat(Zmin,length(Population.objs),1)+1e-10)+1e-10;
    Cosine1   = 1 - pdist2(PopObj1,PopObj1,'cosine');
    Cosine1   = Cosine1.*(1-eye(size(PopObj1,1)));
     Temp1     = sort(Cosine1,2,'descend');
     
    PopObj2 = (ArcPop.objs-repmat(Zmin,length(ArcPop.objs),1))./(repmat(max(AllPop.objs),length(ArcPop.objs),1)-repmat(Zmin,length(ArcPop.objs),1)+1e-10)+1e-10;
    Cosine2   = 1 - pdist2(PopObj2,PopObj2,'cosine');
    Cosine2   = Cosine2.*(1-eye(size(PopObj2,1)));
    Temp2     = sort(Cosine2,2,'descend');
     
    Temp = [Temp1;Temp2];
   % [~,Rank] = sortrows(-Temp);
   % Rank = sort(Rank);
   % [~,Rank] = sort(Temp(:,floor(sqrt(N))+1));
    Rank = Temp(:,floor(sqrt(N))+1);
    
    CV1 = sum(max(0,Population.cons),2);
    CV2 = sum(max(0,ArcPop.cons),2);
    
    Angle1 = Rank(1:N);
    Angle2 = Rank(N+1:length(AllPop));
    
    Index1 = randi(N,1,N);
    Index2 = randi(length(ArcPop),1,N);
    
     i = 0;
    while length(MatingPool)< N  
        
       % if rand < 0.5
      if CV1(Index1(i+1))< CV2(Index2(i+1))     
          MatingPool = [MatingPool,Population(Index1(i+1))];
      else
          MatingPool = [MatingPool,ArcPop(Index2(i+1))];
      end
      
%       if CV1(Index1(i))< CV2(Index2(i))     
%           Matingpool = [MatingPool,Population(Index1(i))];
%       else
%           Matingpool = [MatingPool,ArcPop(Index2(i))];
%       end

      if Angle1(Index1(i+2))<= Angle2(Index2(i+2)) 
          MatingPool = [MatingPool,Population(Index1(i+2))];
      else
          MatingPool = [MatingPool,ArcPop(Index2(i+2))];
      end    
      % else
%       if Angle1(Index1(i+1))< Angle2(Index2(i+1))
%           MatingPool = [MatingPool,Population(Index1(i+1))];
%       else
%           MatingPool = [MatingPool,ArcPop(Index2(i+1))];
%       end 
%        end
      
        i = i + 2 ;
    end
    
    
end

end
