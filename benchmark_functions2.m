
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
    U0 = diag(alpha);

    betta = ones(n,1);
    order=find(sort_idx==count);
    m=size(Y,2);

    order(all(order==0,2),:)=[]
    dim=size(order,1);

    for i=1:dim
        betta(order(i))=x(i);
    end
    Umin = minU0*ones(n, n);

    P = inv(D)*W;
    I= eye(n);
    U = U0+Umin;
    newY = betta.*Y;
    F = pinv(I-P+U)*U*newY;

    [maxF, idF] = max(F,[],2);

    Ug=U0+Umin
    L=D-W
    P=F-betta.*Y
    value = trace(F'*L*F) + trace(P'*Ug*D*P)

    PHI= value;
end
