function Pre=Self_Training(L,t,U,sort_idx,NaNE)

count=1;
Pre=zeros(size(U,1),1);
M_z=[];
while 1
  
    index=find(sort_idx==count);
 
    classifyU=U(index,:);
    KNN_label=KNNC(L,t,classifyU,NaNE);
    Pre(index)=KNN_label;
    for i=1:length(index)
        L=[L;U(index(i),:)];
        t=[t;KNN_label(i)];  
    end
    pos=find(Pre==0);
    M_z(count)=length(pos);
    if length(pos)==0  
        break;
    else  
        if count>2
            if M_z(count-1)==M_z(count)
               sort_idx(pos)=count+1;
            end
        end
    end
    count=count+1;
end
end
