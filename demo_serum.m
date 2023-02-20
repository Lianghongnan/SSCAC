% The code is written by Jie Wen, 
% if you have any problems, please don't hesitate to contact me: wenjie@hrbeu.edu.cn 
% If you find the code is useful, please cite the following reference:
% Jie Wen , Xiaozhao Fang, Yong Xu, Chunwei Tian, Lunke Fei, 
% Low-Rank Representation with Adaptive Graph Regularization [J], 
% Neural Networks, 2018.
% homepage: https://sites.google.com/view/jerry-wen-hit/home
clear all
clc
clear memory;

load('newIdea_spo.mat')
% %name = 'newIdea_spo';
data=spo_fea;
selected_class = 6;
lambda1 = 1e-4;
lambda2 = 1e-2;
lambda3 = 1e-1;
% 
% %load(name)
% fea = double(data);
 nnClass = 6;     % The number of classes
% select_sample = [];
% % select_gnd    = [];
% % for i = 1:selected_class
% %     idx = find(gnd == i);
% %     idx_sample    = fea(idx,:);
% %     select_sample = [select_sample;idx_sample];
% %     %select_gnd    = [select_gnd;gnd(idx)];
% % end
% 
fea = data';
%fea = fea';
fea = fea./repmat(sqrt(sum(fea.^2)),[size(fea,1) 1]);
%gnd = select_gnd;
c   = selected_class;

X = fea;
clear fea select_gnd select_sample idx
[m,n] = size(X);

% ---------- initilization for Z and F -------- %
options = [];
options.NeighborMode = 'KNN';
options.k = 10;
options.WeightMode = 'Binary';      % Binary  HeatKernel
Z = constructW(X',options);
Z = full(Z);
Z1 = Z-diag(diag(Z));         
Z = (Z1+Z1')/2;
DZ= diag(sum(Z));
LZ = DZ - Z;                
[F_ini, ~, evs]=eig1(LZ, c, 0);
Z_ini = Z;
clear LZ DZ Z fea Z1

max_iter= 80;
Ctg = inv(X'*X+2*eye(size(X,2)));

%����õ�FԤ�⣬���ǲ�û��������Ƕ�Z���е�Ncut
[Z,S,U,F,E] = LRR_AGR(X,F_ini,Z_ini,c,lambda1,lambda2,lambda3,max_iter,Ctg);
%% �����õ���Ncut���ѵõ���Z�����׾��ࣩ
addpath('Ncut_9');
% Z2=spo_fea
c=6
Z_out = Z;
A = Z_out;
A = A - diag(diag(A));
A = (A+A')/2;  

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(A,c);
%result_label ��Ncut�� ���
[value,result_label] = max(NcutDiscrete,[],2);
Sil = silhouette(data,result_label);%����ϵ��Ӧ���Ƕ�ԭʼ���ݵģ�����
Sil = mean(Sil)
%result = ClusteringMeasure(gnd, result_label) %result = [ACC MIhat Purity];