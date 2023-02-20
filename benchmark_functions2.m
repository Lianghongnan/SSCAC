
%%
function [PHI]=benchmark_functions2(x,dim,W,Y,sort_idx,fea,count,F)
    [n] = size(W,1);
    maxU0 = 1e5;
    minU0 = 1e-12;
    %%
    D = diag(sum(W));
    alpha = zeros(n,1);
    [label_num,~]=size(find(any(Y,2)));
    alpha (1:label_num,1)=maxU0;
    U0 = diag(alpha);%生成带有未知数的参数矩阵
   %% 插入一个排除正确标签？
    %% 这里也要改成可变的
    betta = ones(n,1);
    order=find(sort_idx==count);
    m=size(Y,2);
%     for i=1:dim
%         if max(F(order(i))) > 1/m
%         order(i)=0;
%         end
%     end
    order(all(order==0,2),:)=[]
    dim=size(order,1);
    %% 以上是新修改的
    for i=1:dim
        betta(order(i))=x(i);
    end
    Umin = minU0*ones(n, n);
    %% 根据这个计算目标函数，得到F值，也就是预测结果
    P = inv(D)*W;
    I= eye(n);
    U = U0+Umin;
    newY = betta.*Y;
    F = pinv(I-P+U)*U*newY;
    %%
%     F = pinv(D+U0-W+Umin)*U0*Y;
    [maxF, idF] = max(F,[],2);
    %%
    Ug=U0+Umin
    L=D-W
    P=F-betta.*Y
    value = trace(F'*L*F) + trace(P'*Ug*D*P)
    %% 以上部分是把本次的轮廓系数计算出来
    PHI= value;
end
