
function [sort_idx,NaNE,NaN,t_index,L,t,U,U_t,dc]=preHandle(fea,gnd,r)
nnClass = length(unique(gnd));  
num_Class=[];
for i=1:nnClass
  num_Class=[num_Class length(find(gnd==i))]; %The number of samples of each class
end
L = [];
t = [];
t_index = [];
 ratio = 0.1;
  for  j=1:nnClass
        idx=find(gnd==j);
        randIdx=randperm(num_Class(j)); 
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
fea2 = fea;
fea2(t_index,:)=[];
U=fea2;
gnd2=gnd;
gnd2(t_index,:)=[];
U_t = gnd2;
data=[L;U];
A=data;
[num,dim]=size(A); 
mdist=pdist(A);
A=tril(ones(num))-eye(num);
[x,y]=find(A~=0);
xx=[x y mdist'];
N=size(xx,1);
percent=2.0;
position=round(N*percent/100);
sda=sort(xx(:,3));
dc=sda(position);
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
k=dc;
[c,~] = size(unique(t));
[arrows,t1,center_idxs]=DPC(data,k,c)
sort_idx = Find_index( arrows,L,U )
[NaNE,NaN]=FindNaN(L,U)
end
