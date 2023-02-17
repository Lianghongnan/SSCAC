function [idF,maxF,valueResult,U0,Y,AccResult,nmiResult]=stdp_lrr_self(Umin,W,D,U0,Y,sort_idx,fea,new_gnd)
%%  Use a self-traning method to predict the class label of unlabeled natural neighbors
%% Initialize parameters
count=1;
maxU0 = 1e5;
samp_num=size(Umin,1);
idF=zeros(size(Umin,1),1);%�����һ��
M_z=[];
labelMatrix=zeros(size(Umin,1),4)
while 1
%%
%     F = inv(D+U0-W+Umin)*U0*Y;
%     F = pinv(D+U0-W+Umin)*U0*Y; 
    %%  С��һ��F�ı�ʾ��ʽ���ѵ���������⣿
    % ���̣�12����Ӧ�Ľ⣨14���ڶ��б�ʾ��ʽ
    P = inv(D)*W;
    I= eye(samp_num);
    U = U0+Umin;
    F = pinv(I-P+U)*U*Y;  
    %% ��ʱ�Ѿ�����˱�ǩ  ��FF��
    [maxF, idF] = max(F,[],2);%idF�������Ҫ��Ԥ����
    for j = 1:samp_num
        FF(j,idF(j)) = 1;%FF����������Y���󣬰�ֵ��ļӽ�ȥ
    end  
    result = ACC_ClusteringMeasure(new_gnd, idF)
    AccResult(count,:)=result;
    z = nmi(new_gnd, idF)
    nmiResult(count,:)=z;
    %% ��¼ÿ�κ���ֵ -->(12)��ֵ
    Ug=U0+Umin;
    L=D-W;
    P=F-Y;
    %%  �ֿ�����
    value1 = trace(F'*L*F);
    value2=trace(P'*Ug*D*P);
    valueResult(count,1)=value1;
    valueResult(count,2)=value2;
    valueResult(count,3)=value1 + value2; %valueResult�ĵ����о���Ҫ��¼�ĺ���ֵ
    %% ��������ó��ޱ������
    pos=find(sort_idx==count);
    if length(pos)==0
        break; %Stopping condition
    end
    index=pos;
    [indexL,~] = size(index);%��ʾ��������Y�ĸ���
%% update U0 and Y
      for i= 1:indexL
            Y(index(i),idF(index(i)))=1;
            labelMatrix(index(i),1:3) = [index(i),new_gnd(index(i)),idF(index(i))];
      end
    %iterLabel = sum(any(Y,2))%�ҵ�Y�в�ȫΪ0������
   labelMatrix(all(labelMatrix==0,2),:)=[]
    [ml,~]=size(labelMatrix)
    %
%%    ����ȥ���Ļ���������ʵ���� 
    ElitistCheck=1; % GSA Parameter--�Ƿ�ʹ����һ�ּ���kbest�ķ���
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