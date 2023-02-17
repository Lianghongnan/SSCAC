function [SSCAC_final]=main(DATA, runtimes)

for r=1:runtimes
gnd=DATA(:,2);
org_fea=DATA(:,3:end);
% org_fea = fea;
[sort_idx,NaNE,NaN,t_index,L,t,U,U_t,dc]=preHandle(org_fea,gnd,r);
new_t = [t;U_t];
%%
new_fea = [L;U];
% reduce demension by PCA
options = [];
options.PCARatio = 0.98;
%options.ReducedDim = 60;
new_fea = knnimpute(new_fea);
[eigvector,eigvalue,meanData,new_data] = PCA(new_fea,options);%将原始数据降维
fea=new_data;
%% 拼接新的索引+标签
label = gnd(t_index)
l_data=1:size(L,1);
l_data = l_data';
label1 = [l_data,label];
%% samp_num 是总的基因数量
samp_num = size(fea,1);
[nnClass,~] = size(unique(label1(:,2)))
num_Class=[];
[nLabel,~] = size(label1); 
fea =  NormalizeFea(fea);
test = fea;
%% 待调节的参数
gama_1 = 100;
gama_2 = 1;
%%
minU0 = 1e-12;
maxU0 = 1e5;
%得到最优解的Z 构建新的权重矩阵W
%%
[A W obj Z] = LRSA(test', gama_1, gama_2);   
A = NormalizeFea(A);%归一化
%     AG = A;
%     save('AG','AG');
%     save('obj','obj');
D= diag(sum(A));
W = A; %这个A就是论文里的W (W=（Z+Z')/2)
%%
%---------------------------------------------------------------   
        Y = zeros(samp_num, nnClass);%samp_num是基因数，nnClass是已知的分为几类
        cLab = zeros(samp_num, nnClass);
        FF = zeros(samp_num, nnClass);
        TestF = ones(samp_num, nnClass);
        U0 = zeros(samp_num, samp_num);
        Umin = minU0*ones(samp_num, samp_num);
 
%%  初始U和Y
    for i = 1:nLabel
        U0(label1(i,1),label1(i,1)) = maxU0;
        Y(label1(i,1),label1(i,2)) = 1;
    end  
    
    %%
%     [ F,Z,W,E] = NNLRR( new_fea',U0,Y,1,145);
    %%
        F = inv(D+U0-W+Umin)*U0*Y;
        [maxF, idF] = max(F,[],2);%idF是最后想要的预测结果
        for j = 1:samp_num
            FF(j,idF(j)) = 1;%FF是重新排列Y矩阵，把值大的加进去
        end
        %% 新的的自训练
        new_gnd=[label;U_t];
        [idF,maxF,valueResult,U0,Y,AccResult,nmiResult]=stdp_lrr_self(Umin,W,D,U0,Y,sort_idx,fea,new_gnd)
        %% 
        %% STDPNF+GFHF
        gfhf_result(r,1)=AccResult(1,1);
        gfhf_result(r,2)=nmiResult(1,1);
        %%
        acc1 = max(AccResult);
        mm = find(AccResult==acc1)
        clear AccResult;
        m(r,1) = mm(1)
        nmi1 = nmiResult(m(r,1),1);
        clear nmiResult;
        result(r,1)=acc1;
        result(r,2)=nmi1;
end
%% 10次平均值

SSCAC_final = mean(result,1)
end