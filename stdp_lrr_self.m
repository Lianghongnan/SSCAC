function [idF,maxF,valueResult,U0,Y,AccResult,nmiResult]=stdp_lrr_self(Umin,W,D,U0,Y,sort_idx,fea,new_gnd)
%%  Use a self-traning method to predict the class label of unlabeled natural neighbors
%% Initialize parameters
count=1;
maxU0 = 1e5;
samp_num=size(Umin,1);
idF=zeros(size(Umin,1),1);%这里改一下
M_z=[];
labelMatrix=zeros(size(Umin,1),4)
while 1
%%
%     F = inv(D+U0-W+Umin)*U0*Y;
%     F = pinv(D+U0-W+Umin)*U0*Y; 
    %%  小修一下F的表示形式，难道是这个问题？
    % 方程（12）对应的解（14）第二行表示形式
    P = inv(D)*W;
    I= eye(samp_num);
    U = U0+Umin;
    F = pinv(I-P+U)*U*Y;  
    %% 此时已经获得了标签  在FF中
    [maxF, idF] = max(F,[],2);%idF是最后想要的预测结果
    for j = 1:samp_num
        FF(j,idF(j)) = 1;%FF是重新排列Y矩阵，把值大的加进去
    end  
    result = ACC_ClusteringMeasure(new_gnd, idF)
    AccResult(count,:)=result;
    z = nmi(new_gnd, idF)
    nmiResult(count,:)=z;
    %% 记录每次函数值 -->(12)的值
    Ug=U0+Umin;
    L=D-W;
    P=F-Y;
    %%  分开计算
    value1 = trace(F'*L*F);
    value2=trace(P'*Ug*D*P);
    valueResult(count,1)=value1;
    valueResult(count,2)=value2;
    valueResult(count,3)=value1 + value2; %valueResult的第三列就是要记录的函数值
    %% 依据序号拿出无标记样本
    pos=find(sort_idx==count);
    if length(pos)==0
        break; %Stopping condition
    end
    index=pos;
    [indexL,~] = size(index);%表示本次重置Y的个数
%% update U0 and Y
      for i= 1:indexL
            Y(index(i),idF(index(i)))=1;
            labelMatrix(index(i),1:3) = [index(i),new_gnd(index(i)),idF(index(i))];
      end
    %iterLabel = sum(any(Y,2))%找到Y中不全为0的行数
   labelMatrix(all(labelMatrix==0,2),:)=[]
    [ml,~]=size(labelMatrix)
    %
%%    这里去掉的话就是消融实验了 
    ElitistCheck=1; % GSA Parameter--是否使用另一种计算kbest的方法
    [Fbest,Lbest,BestChart]=GSA(indexL,ElitistCheck,W,Y,sort_idx,fea,count,F)
%
    for i= 1:indexL
        Y(index(i),idF(index(i)))=Lbest(i);
        labelMatrix(i,4)=Lbest(i);
    end
    
   save(strcat('D:\MATLAB\R2019b\bin\LRR-Filter-GA-self\Confidence_reesult\',num2str(count)),'labelMatrix');
   for i = 1:ml     
         U0(labelMatrix(i,1),labelMatrix(i,1)) = maxU0;
   end
%%
  clear labelMatrix;
%% 
 count=count+1;
end
end