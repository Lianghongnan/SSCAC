
function [X]=initialization(dim,N,up,down)
if size(up,1)>1
    for i=1:dim
    high=up(i);low=down(i);
    X(i,:)=rand(1,N).*(high-low)+low;
    end
end
s=up*dim;
if size(up,1)==1
    for i = 1:N
        X(:,i)= (makeGDSum(s,dim))
    end
end
