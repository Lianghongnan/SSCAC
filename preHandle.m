
%% 数据预处理，把这个和整体整合在一起
function [sort_idx,NaNE,NaN,t_index,L,t,U,U_t,dc]=preHandle(fea,gnd,r)
%%
load ( 'D:\MATLAB\R2019b\bin\LRR-Filter-GA-self\有标签数据（非基因）\seed_value\seed_ecoli.mat');
nnClass = length(unique(gnd));  % The number of classes;
num_Class=[];
for i=1:nnClass
  num_Class=[num_Class length(find(gnd==i))]; %The number of samples of each class
end

%%
L = [];
t = [];
t_index = [];
 %% 这里写一个按照类标签比例选取的
 ratio = 0.1;%每次更改是改这里
  for  j=1:nnClass
        idx=find(gnd==j);
        % randperm(n) 返回行向量，其中包含从 1 到 n 没有重复元素的整数随机排列。
        rand('seed',3);
        %% 随机数
        %rand('seed',seed_e(r));%固定随机种子
%         s=27+r;
%         rand('seed',s);%固定随机种子
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
%% 获取LRRADP需要的数据
data=[L;U];
% 想在这里直接计算截断距离
A=data;
[num,dim]=size(A); 
mdist=pdist(A);%两两行之间距离
A=tril(ones(num))-eye(num);
[x,y]=find(A~=0);
xx=[x y mdist'];
N=size(xx,1);
percent=2.0;
position=round(N*percent/100);
sda=sort(xx(:,3));%所有距离的升序排列
dc=sda(position);%截断距离选取的所有距离的2%
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
%%
k=dc;
% c=4;%分成几类
[c,~] = size(unique(t));%自动确定共有几类
[arrows,t1,center_idxs]=DPC(data,k,c)
sort_idx = Find_index( arrows,L,U )
[NaNE,NaN]=FindNaN(L,U)
end
