function [idF,maxF,valueResult,U0,Y,AccResult,nmiResult]=stdp_lrr_self(Umin,W,D,U0,Y,sort_idx,fea,new_gnd)
count=1;
maxU0 = 1e5;
samp_num=size(Umin,1);
idF=zeros(size(Umin,1),1);
M_z=[];
labelMatrix=zeros(size(Umin,1),4)
while 1
    P = inv(D)*W;
    I= eye(samp_num);
    U = U0+Umin;
    F = pinv(I-P+U)*U*Y;  
    [maxF, idF] = max(F,[],2);
    for j = 1:samp_num
        FF(j,idF(j)) = 1;
    end  
    result = ACC_ClusteringMeasure(new_gnd, idF)
    AccResult(count,:)=result;
    z = nmi(new_gnd, idF)
    nmiResult(count,:)=z;
    Ug=U0+Umin;
    L=D-W;
    P=F-Y;

    value1 = trace(F'*L*F);
    value2=trace(P'*Ug*D*P);
    valueResult(count,1)=value1;
    valueResult(count,2)=value2;
    valueResult(count,3)=value1 + value2;
    pos=find(sort_idx==count);
    if length(pos)==0
        break; 
    end
    index=pos;
    [indexL,~] = size(index);
      for i= 1:indexL
            Y(index(i),idF(index(i)))=1;
            labelMatrix(index(i),1:3) = [index(i),new_gnd(index(i)),idF(index(i))];
      end
   labelMatrix(all(labelMatrix==0,2),:)=[]
    [ml,~]=size(labelMatrix)

    ElitistCheck=1;
    [Fbest,Lbest,BestChart]=GSA(indexL,ElitistCheck,W,Y,sort_idx,fea,count,F)

    for i= 1:indexL
        Y(index(i),idF(index(i)))=Lbest(i);
        labelMatrix(i,4)=Lbest(i);
    end
   for i = 1:ml     
         U0(labelMatrix(i,1),labelMatrix(i,1)) = maxU0;
   end

  clear labelMatrix;
 count=count+1;
end
end
