
%This function initializes the position of the agents in the search space, randomly.
function [X]=initialization(dim,N,up,down)

% if size(up,1)==1
%     X=rand(dim,N).*(up-down)+down;
% end
%% ��������dim*N
% ����ÿ���ڼ����ʼ�������ֶ�����������

if size(up,1)>1
    for i=1:dim
    high=up(i);low=down(i);%����high=1,low = 0.09
    X(i,:)=rand(1,N).*(high-low)+low;
    end
end
s=up*dim;
if size(up,1)==1
    for i = 1:N
        %% makeGDSum�ý����һ�����֣���Ҫ�����У�������ҪN��
        X(:,i)= (makeGDSum(s,dim))
    end
end