clc,clear;
folder_now = pwd;
addpath([folder_now, '\funs']);
addpath 'D:\MATLAB\R2019b\bin\Code for NNLRR\YALL1_v1.3'
runtimes = 100;
load ('.\data\gal.mat', 'new_Gal');

DATA=new_Gal;
for r=1:runtimes
gnd=DATA(:,2);
org_fea=DATA(:,3:end);
[sort_idx,NaNE,NaN,t_index,L,t,U,U_t,dc]=preHandle(org_fea,gnd,r);
new_t = [t;U_t];
%%
new_fea = [L;U];
options = [];
options.PCARatio = 0.98;
new_fea = knnimpute(new_fea);
[eigvector,eigvalue,meanData,new_data] = PCA(new_fea,options);
fea=new_data;
label = gnd(t_index)
l_data=1:size(L,1);
l_data = l_data';
label1 = [l_data,label];
samp_num = size(fea,1);
[nnClass,~] = size(unique(label1(:,2)))
num_Class=[];
[nLabel,~] = size(label1); 
fea =  NormalizeFea(fea);
test = fea;
gama_1 = 100;
gama_2 = 1;
minU0 = 1e-12;
maxU0 = 1e5;
[A W obj Z] = LRSA(test', gama_1, gama_2);   
A = NormalizeFea(A);
D= diag(sum(A));
W = A; 
        Y = zeros(samp_num, nnClass);
        cLab = zeros(samp_num, nnClass);
        FF = zeros(samp_num, nnClass);
        TestF = ones(samp_num, nnClass);
        U0 = zeros(samp_num, samp_num);
        Umin = minU0*ones(samp_num, samp_num);
    for i = 1:nLabel
        U0(label1(i,1),label1(i,1)) = maxU0;
        Y(label1(i,1),label1(i,2)) = 1;
    end  
        F = inv(D+U0-W+Umin)*U0*Y;
        [maxF, idF] = max(F,[],2);
        for j = 1:samp_num
            FF(j,idF(j)) = 1;
        end
        new_gnd=[label;U_t];
        [idF,maxF,valueResult,U0,Y,AccResult,nmiResult]=stdp_lrr_self(Umin,W,D,U0,Y,sort_idx,fea,new_gnd)
        acc1 = max(AccResult);
        mm = find(AccResult==acc1)
        clear AccResult;
        m(r,1) = mm(1)
        nmi1 = nmiResult(m(r,1),1);
        clear nmiResult;
        result(r,1)=acc1;
        result(r,2)=nmi1;
end
SSCAC_final = mean(result,1)


