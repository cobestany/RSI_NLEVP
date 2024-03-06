function [A,M] = genAM(Prob)
%% [A,M] = genAM(Prob) returns matrices A and M
%   of the linearized eigenvalue problem:
%               A*w = lam*M*w
%   where lam has the same spectrum as the original NLEVP.
%

%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials
n = size(Ci{1},1);  %--- problem size

%--- construct A
I = eye(n);
%---------- last row of A
%---------- left part corresponding to rational approx
B = [];
for ik = 1:np
    B = [B Bi{ik}];
end
%---------- right part corresponding to Chebyshev approx
%---------- scale
scal = 2*c/r;
%---------- special cases
if deg == 0
    C = Ci{1};
    Atop = [diag(zk) ones(np,1)];
    Abot = [];
elseif deg == 1
    C = Ci{1} - c/r*Ci{2};
    Atop = [diag(zk) ones(np,1)];
    Abot = [];
else
    C = [];
    for ik = 1:deg-2
        C = [C Ci{ik}];
    end
    C = [C Ci{deg-1}-Ci{deg+1} Ci{deg}-scal*Ci{deg+1}];
    %---------- top part corresponding to rational approx
    Atop = [diag(zk) ones(np,1) zeros(np,deg-1)];
    %---------- bottom part corresponding to Chebyshev approx
    Abot = [zeros(deg,np) diag(ones(deg-1,1),-1)+scal*diag(ones(deg,1))+diag(ones(deg-1,1),1)];
    Abot(1,np+1) = Abot(1,np+1)/2;
    Abot = Abot(1:end-1,:);
end
A = [Atop; Abot];
A = kron(A,I);
A = [A;B C];

%--- construct M
%---------- special cases
if deg == 0
    t0 = [ones(1,np), 0];
    M = kron(diag(t0),I);
elseif deg == 1
    t0 = [ones(1,np), 0];
    t1 = [zeros(1,np), -1/r];
    M = kron(diag(t0),I)+kron(diag(t1),Ci{deg+1});
else
    t0 = [ones(1,np), 1/r, zeros(1,deg-1)];
    t2 = [zeros(1,np+deg-1), -2/r];
    t1 = ~(t0+t2);
    M  = kron(diag(t0),I)+kron(diag(t1),(2/r)*I)+kron(diag(t2),Ci{deg+1});
end