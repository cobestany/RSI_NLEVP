function [v] = ShiInv(Prob,mu,v,ninv,nref,f,LUPQ0,LUPQc)
%% [v] = ShiInv(Prob,mu,v,ninv,nref,f,LUPQ0,LUPQc) 
%   computes shift and invert using a given shift 
%   mu and vector v.
% INPUT:
%   mu      - a given shift
%   v       - starting vector
%   ninv    - maximum num. of iterations
%   nref    - num. of iterative refinement
%   f       - 0 for using pre-computed LU;
%             1 for using new LU w.r.t. mu
%   LUPQ0   - pre-computed LU of Schur complement
%   LUPQc   - LU of Ci{end}
% NOTE: Using pre-computed LU factorization unless
%   the shift is far away from the cluster center. 
%   The result can be refined via iterative refinement 
%   and Rayleigh quotient.
%%  iterative refinement x{0} = w
%------ step 1: r{m} = M*w - (A - mu{m}*M)*x{m} = M*(w - M\(A - mu{m}*M)*x{m})
%------ step 2: solve (A - mu{0}*M)*dx{m} = r{m}
%------ step 3: x{m+1} = x{m} + dx{m}


%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials
n = size(Ci{1},1);  %--- problem size

%------------- generate long vector w from v
mu0 = RayQuo(Prob,v);
w = LongVec(Prob,mu0,v);

%--- shift-and-invert
for iter = 1:ninv
    w0 = w;
    v0 = v;
    
    %--- compute (A-mu*M)\(M*w)
    [w,LUPQ0] = ShiInv_solve(Prob,mu,w,LUPQ0,f);
    v = w(np*n+1:np*n+n);
    %--- do not change mu, use the same LU
    f = false;
    
    %--- iterative refinement
    w1 = w;
    v1 = v;
    for ref = 1:nref
        mu1 = RayQuo(Prob,v1);
        w_tmp = ShiInv_solve_rev(Prob,mu1,w1,LUPQc);
        [dw,~] = ShiInv_solve(Prob,mu,w1-w_tmp,LUPQ0,false);
        w1 = w1 + dw;
        v1 = w1(np*n+1:np*n+n);
    end
    v = v1;
    
    %--- safe-guard
    if 1-abs(v0'*v)/norm(v0)/norm(v) < 1e-5   %---- converged
        break
    end
    
    % %--- Rayleigh quotient refinement
    % if mod(iter,fix(ninv/4)) == 0   %---- refine mu and w
    %     lam = RayQuo(Prob,v);
    %     w = LongVec(Prob,lam,v);
    % end
    
    %--- generate w for next iterate
    w = w/norm(v);
end
% fprintf(repmat('\b',1,iter));
fprintf(' ')