function [A W obj Z] = LRSA(X, g1, g2)
% y: num*1 clbetaHter indicator vector
% A: num*num learned Hymmetric Himilarity matrix
addpath(genpath('.\YALL1_v1.3'));
NITER = 1000;
dim = size(X,1);
num = size(X,2);
if nargin < 3
 g2 = 0.1;
end
if nargin < 2
 g1 = 0.1;
end
distX = L2_distance_1(X,X);
islocal = 0;
if islocal == 1
    k = 50;
    [distX1, idx] = sort(distX,2);
    A = 10000000*ones(num);
    for i = 1:num
        di = distX1(i,2:k+2);
        id = idx(i,2:k+2);
        A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
    end;
    A0 = (A+A')/2;
else
    A0 = distX;
end
Z = zeros(num);
E = sparse(dim, num);
H = Z;
gama_1 = g1;
gama_2 = g2;
beta = 25;
Y1 = zeros(dim, num);
Y2 = zeros(num);
normX = norm(X,2);
Ita_a = 2*normX;
ro = 1.1;
svp = 5;
tol = 1e-6;
max_beta = 1e5;
for iter = 1:NITER

    Z_old = Z;
    E_old = E;
    
    Q1 = Z + (X'*(X-X*Z-E+Y1/beta) - (Z - H + Y2/beta))/Ita_a;
    [U D V] = svd(Q1,'econ');
    D = diag(D);
    svp = length(find(D>1/(Ita_a*beta)));
    if svp >= 1
        D = D(1:svp)-1/(Ita_a*beta);
    else
        D = 0;
        svp = 1;
    end
    Z = U(:,1:svp)*diag(D)*V(:,1:svp)';
    
    E1 = X-X*Z+Y1/beta;
    theta = gama_1/beta;
    E = max(0,E1-theta)+min(0,E1+theta);

    H = Calc_H(A0, Z, Y2, gama_2, beta);
    errH = norm(Z-Z_old,2)/normX;
    errE = norm(E-E_old,2)/normX;
    convergence = (errE < tol)&&(errH < tol);
    if convergence
        break;
    end    
    obj(iter) = rank(Z) + gama_1*sum(sum(abs(E))) + gama_2*sum(sum(distX.*Z));
    Y1 = Y1 + beta*(X-X*Z-E);
    Y2 = Y2 + beta*(Z-H);
    beta = min(max_beta, ro*beta);
end

W=(Z+Z')/2;
A = (H+H')/2;


function [A] = Calc_H(disX,Z,Y,gama, miu)
ss = size(disX, 1);
T = 0.0005;
for i=1:ss
    idxa0 = 1:ss;
    
    dz = Z(i,idxa0);
    dy = Y(i,idxa0);
    dx = disX(i, idxa0);
    ad = dz+dy/miu-gama*dx/miu;
    A(i,idxa0) = EProjSimplex_new(ad);
end;



