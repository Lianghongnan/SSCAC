
%% ����Ԥ���������������������һ��
function [sort_idx,NaNE,NaN,t_index,L,t,U,U_t,dc]=preHandle(fea,gnd,r)
%%
load ( 'D:\MATLAB\R2019b\bin\LRR-Filter-GA-self\�б�ǩ���ݣ��ǻ���\seed_value\seed_ecoli.mat');
nnClass = length(unique(gnd));  % The number of classes;
num_Class=[];
for i=1:nnClass
  num_Class=[num_Class length(find(gnd==i))]; %The number of samples of each class
end

%%
L = [];
t = [];
t_index = [];
 %% ����дһ���������ǩ����ѡȡ��
 ratio = 0.1;%ÿ�θ����Ǹ�����
  for  j=1:nnClass
        idx=find(gnd==j);
        % randperm(n) ���������������а����� 1 �� n û���ظ�Ԫ�ص�����������С�
        rand('seed',3);
        %% �����
        %rand('seed',seed_e(r));%�̶��������
%         s=27+r;
%         rand('seed',s);%�̶��������
        randIdx=randperm(num_Class(j)); %randIdx create m random number, m is the size of idx.
        sele = round(num_Class(j) * ratio);
        if num_Class(j) < sele 
            L=[L;fea(idx(randIdx),:)];
            t=[t;gnd(idx(randIdx))];
            t_index = [t_index;idx(randIdx)];
        else
            L=[L;fea(idx(randIdx(1:sele)),:)];
            t=[t;gnd(idx(randIdx(1:sele)))];
            t_index = [t_index;idx(randIdx(1:sele))];
        end
 end
%% 
fea2 = fea;
fea2(t_index,:)=[];
U=fea2;
%%
gnd2=gnd;
gnd2(t_index,:)=[];
U_t = gnd2;
%% ��ȡLRRADP��Ҫ������
data=[L;U];
% ��������ֱ�Ӽ���ضϾ���
A=data;
[num,dim]=size(A); 
mdist=pdist(A);%������֮�����
A=tril(ones(num))-eye(num);
[x,y]=find(A~=0);
xx=[x y mdist'];
N=size(xx,1);
percent=2.0;
position=round(N*percent/100);
sda=sort(xx(:,3));%���о������������
dc=sda(position);%�ضϾ���ѡȡ�����о����2%
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
%%
k=dc;
% c=4;%�ֳɼ���
[c,~] = size(unique(t));%�Զ�ȷ�����м���
[arrows,t1,center_idxs]=DPC(data,k,c)
sort_idx = Find_index( arrows,L,U )
[NaNE,NaN]=FindNaN(L,U)
end
