
%This function initializes the position of the agents in the search space, randomly.
function [X]=initialization(dim,N,up,down)

% if size(up,1)==1
%     X=rand(dim,N).*(up-down)+down;
% end
%% 这个结果是dim*N
% 或者每次在计算初始化的是手动调节上下限

if size(up,1)>1
    for i=1:dim
    high=up(i);low=down(i);%这里high=1,low = 0.09
    X(i,:)=rand(1,N).*(high-low)+low;
    end
end
s=up*dim;
if size(up,1)==1
    for i = 1:N
        %% makeGDSum得结果是一行数字，需要换成列，并且我要N列
        X(:,i)= (makeGDSum(s,dim))
    end
end