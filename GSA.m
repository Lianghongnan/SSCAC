% Gravitational Search Algorithm.
function [Fbest,Lbest,BestChart]=GSA(dim,ElitistCheck,W,Y,sort_idx,fea,count,F)
 Rpower=1;    
 Rnorm=2; 
 N = 50; 
down=ones(dim,1)*0.009;         % Lower Bound of Variables
up= ones(dim,1)*1;
X=initialization(dim,N,up,down); 
K=round(dim/2);
BestChart=[];

V=zeros(dim,N);
Max_Iteration  = 100;
min_flag=1;
for iteration=1:Max_Iteration

    X=space_bound(X,up,down); 
    
    %Evaluation of agents. 
    fitness= evaluateF(X,dim,W,Y,sort_idx,fea,count,F)
   
    if min_flag==1
    [best, best_X]=min(fitness); 
    [best best_X]=max(fitness); 
    end        
    
    if iteration==1
       Fbest=best;Lbest=X(:,best_X);
    end
    if min_flag==1
      if best<Fbest 
       Fbest=best;Lbest=X(:,best_X);
      end
    else 
      if best>Fbest  
       Fbest=best;Lbest=X(:,best_X);
      end
    end
      
BestChart=[BestChart Fbest];

[M]=massCalculation(fitness,min_flag); 

G=Gconstant(iteration,Max_Iteration); 
a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,Max_Iteration);
[X,V]=move(X,a,V);

end 
