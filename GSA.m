% GSA code v1.1.
% Generated by Esmat Rashedi, 2010. 
% "	E. Rashedi, H. Nezamabadi-pour and S. Saryazdi,
% 揋SA: A Gravitational Search Algorithm?, Information sciences, vol. 179,
% no. 13, pp. 2232-2248, 2009."

% Gravitational Search Algorithm.
function [Fbest,Lbest,BestChart]=GSA(dim,ElitistCheck,W,Y,sort_idx,fea,count,F)

%V:   Velocity.
%a:   Acceleration.
%M:   Mass.  Ma=Mp=Mi=M;
%dim: Dimension of the test function.
%N:   Number of agents.
%X:   Position of agents. dim-by-N matrix.
%R:   Distance between agents in search space.
%[low-up]: Allowable range for search space.
%Rnorm:  Norm in eq.8.
%Rpower: Power of R in eq.7.
 Rpower=1;       % GSA Parameter---公式中的一个平方数
 Rnorm=2; %欧式距离的范数
 N = 50; 
%get allowable range and dimension of the test function.

% %%
down=ones(dim,1)*0.009;         % Lower Bound of Variables
up= ones(dim,1)*1;
%%
% down = 1;
% up =adjustU0
%random initialization for agents.
X=initialization(dim,N,up,down); 
K=round(dim/2);
% K=round(dim*0.75)
%create the best so far chart and average fitnesses chart.
BestChart=[];

V=zeros(dim,N);
Max_Iteration  = 100;% 画图的横轴
min_flag=1;
for iteration=1:Max_Iteration
%     iteration
    
    %Checking allowable range. 
    X=space_bound(X,up,down); 
    
    %Evaluation of agents. 
    fitness= evaluateF(X,dim,W,Y,sort_idx,fea,count,F)
    %%
    if min_flag==1
    [best, best_X]=min(fitness); %minimization.
    else
    [best best_X]=max(fitness); %maximization.
    end        
    
    if iteration==1
       Fbest=best;Lbest=X(:,best_X);
    end
    if min_flag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(:,best_X);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(:,best_X);
      end
    end
      
BestChart=[BestChart Fbest];

%Calculation of M. eq.14-20
% 根据适应度计算惯性质量
[M]=massCalculation(fitness,min_flag); 

%Calculation of Gravitational constant. eq.13.
G=Gconstant(iteration,Max_Iteration); 

%Calculation of accelaration in gravitational field. eq.7-10,21.
a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,Max_Iteration);

%Agent movement. eq.11-12
[X,V]=move(X,a,V);

end %iteration

