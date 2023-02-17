
function  fitness=evaluateF(X,dim,W,Y,sort_idx,fea,count,F)
[dim,N]=size(X);
for i=1:N 
%     L=X(:,i)';
    %calculation of objective function
%     [PHI]=benchmark_functions(X,dim,W,Y,sort_idx,fea,count);
    [PHI]=benchmark_functions2(X,dim,W,Y,sort_idx,fea,count,F)
    fitness(i)=PHI;
end