function [res_vec,res_nrm] = get_res(Prob,Lam,U)
%% [res_vec,res_nrm] = get_res(Prob,Lam,U) returns
%   the residual vectors and norm of given eigenpairs:
%               res_vec = T(lam)*U
%               res_norm = norm(res_vec)
%   Here we consider the approximated NLEVP. It can be
%   replaced by the original NLEVP if necessary.
% INPUT:
%   Lam - set of eigenvalues
%   U   - set of eigenvectors (same order as lam)
% Note:
%   NO scaling on res_nrm and NO normalization on res_vec.

%--- init
res_vec = [];
res_nrm = [];
%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials

%--- compute residual
for k = 1:length(Lam)
    %--- load eigenpair
    lam = Lam(k);
    u = U(:,k)/norm(U(:,k));
    res = 0;
    %--- rational part
    for j = 1:np
        res = res + Bi{j}*u/(lam - zk(j));
    end
    %--- Chebyshev part
    if deg == 0
        res = res + Ci{1}*u;
    else
        t0 = u;
        t1 = (lam-c)/r*u;
        res = res + Ci{1}*t0 + Ci{2}*t1;
        for j = 2:deg
            t2 = 2*(lam-c)/r*t1 - t0;
            t0 = t1;
            t1 = t2;
            res = res + Ci{j+1}*t2;
        end
    end
    %--- write residual
    res_vec(:,k) = res;
    res_nrm(k) = norm(res);
end