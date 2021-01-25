function MatingPool = MatingSelection(Population,ArcPop,N)
% SelectedIndex = zeros(1,N);
  MatingPool = [];
if length(ArcPop)<N
    SelectedIndex = TournamentSelection(2,N,-sum(max(0,Population.cons),2));
    MatingPool = Population(SelectedIndex);
else
    AllPop = [ Population,ArcPop]; 
    Zmin       = min(AllPop.objs,[],1);
   % Zmax = max(Population.objs,[],1);
    PopObj = (AllPop.objs-repmat(Zmin,length(AllPop.objs),1))./(repmat(max(AllPop.objs),length(AllPop.objs),1)-repmat(Zmin,length(AllPop.objs),1)+1e-10)+1e-10;
    Cosine   = 1 - pdist2(PopObj,PopObj,'cosine');
    Cosine   = Cosine.*(1-eye(size(PopObj,1)));
    
     Temp     = sort(-Cosine,2);
    [~,Rank] = sortrows(Temp);
    
   
    
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

      if Angle1(Index1(i+2))< Angle2(Index2(i+2))
          MatingPool = [MatingPool,Population(Index1(i+2))];
      else
          MatingPool = [MatingPool,ArcPop(Index2(i+2))];
      end    
     %   else
%       if Angle1(Index1(i+1))< Angle2(Index2(i+1))
%           MatingPool = [MatingPool,Population(Index1(i+1))];
%       else
%           MatingPool = [MatingPool,ArcPop(Index2(i+1))];
%       end 
      %  end
      
        i = i + 2 ;
    end
    
    
end

end
