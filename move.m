%This function updates the velocity and position of agents.
function [X,V]=move(X,a,V)
%movement.
[dim,N]=size(X);
V=rand(dim,N).*V+a; %eq. 11.---�ٶ�---�������������ôд�ģ��ҿ�����
X=X+V; %eq. 12.---����
%%
% out_abs=rand(1,N);
% mid_abs=sum(out_abs)*S_abs;
% out_abs=out_abs/sum(out_abs)*S_abs;
%%  �����Ƕ�ƽ�����Լ��
% S_abs=K;
% for i = 1: N
%     mid_x(:,i) = sum(X(:,i));
%     final_X(:,i)=X(:,i)/mid_x(:,i)*S_abs
% end
% 
% end