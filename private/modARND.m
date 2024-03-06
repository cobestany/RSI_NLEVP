function [lam,U] = modARND(Prob,mu,nit)
%% Modified Arnoldi method
%   Arnoldi method trimmed for large generalized EVP:
%                A*w = lam*M*w
%   Here we consider either storing short vectors ONLY.
% INPUT:
%   mu    - shift, cannot have mu == c
%   nit   - num. of iterates
    
%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials
n = size(Ci{1},1);  %--- problem size

%% ------------ short vec --------------
%%---- init
v0 = randn(n,1);
v0 = v0/norm(v0);
V = v0;
H = zeros(nit+1,nit);
if abs(mu - c) < 1e-12
    mu = c + r*1e-5*randn(1,'like',1i);
end
LUPQ0.L = []; LUPQ0.U = []; LUPQ0.P = []; LUPQ0.Q = [];
f = 1;
%%---- ARND iteration
for k = 2:nit+1
    %------- generate long vectors
    v0 = V(:,k-1);
    mu0 = RayQuo(Prob,v0);
    q0 = LongVec(Prob,mu0,v0);
    %------- compute q1 <- M\A*q0
%     q1 = solveLinInv(q0,mu,Bi,Ci,c,r,zk,Lc,Uc);
    %------- compute q1 <- (A-mu*M)\M*q0
    [q1,LUPQ0] = ShiInv_solve(Prob,mu,q0,LUPQ0,f);
    f = 0;
    v1 = q1(np*n+1:np*n+n);
    %------- modified GS
    for j = 1:k-1
        vj = V(:,j);
        t = vj'*v1;
        H(j,k-1) = t;
        v1 = v1 - t*vj;
    end
    %------- normalization
    t = norm(v1,2);
    H(k,k-1) = t;
    if (t < 1e-10),break,end
    t = 1.0/t;
    V = [V, v1*t];
    fprintf('-')
end

%% ------------ solve --------------
nn = min(nit,size(V,2));
V = V(:,1:nn);
Hm = H(1:nn,1:nn);
[X,eval] = eig(Hm);
% lam = diag(eval);
lam = 1./diag(eval)+mu;
U = V*X;
%------------ exclude points outside of region of interest
idx = abs(lam-c)<0.8*r;
lam = lam(idx);
U = U(:,idx);
fprintf('\n')