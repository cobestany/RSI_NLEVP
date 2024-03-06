function [eval,evec,Prob,shift] = RSI(T,params,nrm)
%% [eval,evec,Prob] = RSI(T,params,nrm) 
% Reduced Subspace Iteration:
%   for a matrix-valued nonlinear function T, parameters params, and 
%   (optional) a scaling function nrm, returns the eigenvalues eval and
%   eigenvectors evec of the nonlinear eigenvalue problem T(lam)*u = 0.
%
% INPUTS:
%   T      -> a matrix-valued nonlinear function handler, T(z) is a matrix
%   params -> a struct of parameters
%     * params.deg       -> degree of Chebyshev polynomials
%     * params.n_quad    -> number of quadrature points
%     * params.center    -> center of region of interest
%     * params.radius    -> radius of region of interest
%     * params.tol       -> tolerance of relative residual
%     * params.nit_arnd  -> number of iterations for surrogate Arnoldi
%                           (default: 30)
%     * params.n_shift   -> number of initial shifts (default: 1)
%     * params.n_restart -> number of restart per subregion (default: 1)
%     * params.nit_inv   -> number of iterations for shift-and-invert
%                           (default: 5)
%     * params.nit_ref   -> number of iterations for iterative refinement
%                           (default: 0)
%     * params.dim       -> initial dimension of projection subspace
%                           (default: 10)
%     * params.buffer    -> split subregion if get more than (dim - buffer) 
%                           eigenvaluess in one subregion (default: dim/2)
%   nrm    -> (optional) a scalar function handler for scaling relative 
%             residual, e.g., rel_res = norm(T(lam)*u)/norm(u)/nrm(lam).
%             If not given, nrm(z) = 1.0 by default.
%
% OUTPUTS:
%   eval   -> eigenvalues stored in a vector
%   evec   -> corresponding eigenvectors stored in columns of evec
%   Prob   -> (optional) a struct of Chebyshev-rational approximation info
%   shift  -> (optional) all shifts for shift-and-invert
%
% NOTE:
%   Currently the region of approximation and the region of interest are 
%   both circles with the same center. However, the region of interest can 
%   be any shape fully enclosed by the region of approximation.

%% Extract parameters
deg = params.deg;     %--- degree of Chebyshev polynomials
np  = params.n_quad;  %--- number of quadrature points
c   = params.center;  %--- center of region of approximation
r   = 1.25*params.radius;  %--- radius of region of approximation
tol = params.tol;     %--- tolerance of relative residual
if ~isfield(params,'nit_arnd')
    nit = 30;
else
    nit = params.nit_arnd;
end
if ~isfield(params, 'n_shift')
    m = 1;
else
    m = params.n_shift;
end
if ~isfield(params, 'n_restart')
    nrst = 1;
else
    nrst = params.n_restart;
end
if ~isfield(params, 'nit_inv')
    ninv = 5;
else
    ninv = params.nit_inv;
end
if ~isfield(params, 'nit_ref')
    nref = 0;
else
    nref = params.nit_ref;
end
if ~isfield(params, 'dim')
    ns_init = 10;
else
    ns_init = params.dim;
end
if ~isfield(params, 'buffer')
    nsp = fix(ns_init/2);
else
    nsp = params.buffer;
end
%--- scaling func for relative residual
if nargin == 2, nrm = @(z) 1.0;  end

%% Chebyshev-Rational Approximation
%---------- get Chebyshev part coefficient matrices
[Ci] = cheb_approx(T,deg,c,r);

%---------- get rational part coefficient matrices
[zk, omega] = contQuad(np,c,1.25*r,0);   %--- quadrature points and weights
Bi = {};
for k = 1:np
    z = zk(k);
    Tz = T(z);
    Tcheb = cheb_val(Ci,c,r,z);
    Bi{k} = omega(k)*(Tz-Tcheb);
end

%% Initialization
%---------- construct problem
Prob = struct('Bi',{Bi},'Ci',{Ci},'c',c,'r',r,'sig',zk,'np',np,'deg',deg);

%---------- initial LU factorizations
[LUPQc.L,LUPQc.U,LUPQc.P,LUPQc.Q] = lu(Ci{end});

%---------- compute shifts
[z,~] = modARND(Prob,c,nit);
[~,shift]= Kmeans(m,z);  %--- initial shifts
shift_idx = 1;
current_shift_idx = 0;

%---------- merge Clusters if two are close
merge_cluster = true;
while merge_cluster
    [d,i,j] = minDist(shift);
    if d < 0.1*r
        fprintf('merge subregions...\n')
        new_shift = (shift(i)+shift(j))/2;
        shift([i j]) = [];
        shift = [shift; new_shift];
        m = m-1;
    else
        merge_cluster = false;
    end     
end

%% Reduced Subspace Iteration
eval = [];      %--- locked eigenvalues
evec = [];      %--- locked eigenvectors
ev   = [];      %--- new eigenvalues
V    = [];      %--- basis vectors
Vs   = [];      %--- Schur vectors of locsked eigenvectors
rst  = 1;       %--- counter for restarts
ns   = ns_init; %--- initial dimension of projection subspace
while shift_idx <= length(shift)
    fprintf('cluster: %d/%d\n',shift_idx, length(shift))
    %---------- generate basis V
    l = length(ev);
    for k = 1:ns
        if current_shift_idx ~= shift_idx
            current_shift_idx = shift_idx;
            mu0 = shift(shift_idx);
            LUPQ0 = initLU(Prob,mu0);
        end
        if k <= l
            v = U(:,k);
            mu = ev(k);
        else
            v = randn(size(Ci{1},2),1);
            mu = mu0;
        end
        if abs(mu-mu0)<0.01*r
            f = false;
        else
            f = true;
        end
        v = ShiInv(Prob,mu,v,ninv,nref,f,LUPQ0,LUPQc);  
        V = mgs(V,k,v);
        fprintf('(%d/%d)\n',k,ns);
    end
    %---------- solve projected NLEVP
    Bip = {};
    for k = 1:np
        Bip{k} = V'*(Bi{k}*V);
    end
    Cip = {};
    for k = 1:deg+1
        Cip{k} = V'*(Ci{k}*V);
    end
    projProb = struct('Bi',{Bip},'Ci',{Cip},'c',c,'r',r,'sig',zk,...
                    'np',np,'deg',deg);
    [ev,U] = solveNLEVP(projProb);
    [ev,idx] = in_cluster(ev,shift,shift_idx);
    U = V*U(:,idx);
    %---------- compute residual and lock eigenpairs
    res_nrm = [];
    for i = 1:length(ev)
        u = U(:,i);
        res_vec = T(ev(i))*u/norm(u);
        res_nrm(i) = norm(res_vec)/nrm(ev(i));
    end
    idx = find(res_nrm < tol);
    for i = 1:length(idx)
        ui = U(:,idx(i));
        d = 1.0; %--- cosine distance to locked Schur vectors
        if ~isempty(Vs)
            d = norm(ui - Vs*(Vs'*ui))/norm(ui);
        end
        if d < 1e-3,continue,end
        Vs = mgs(Vs,size(Vs,2)+1,ui);
        eval = [eval; ev(idx(i))];
        evec = [evec, ui];
    end
    %---------- split cluster if too many eval
    if length(ev) > ns_init-nsp
        [~,new_shift]= Kmeans(2,ev);
        new_shift = new_shift(~isnan(new_shift));
        shift(shift_idx) = new_shift(1);
        if length(new_shift) > 1
            shift = [shift; new_shift(2:end)];
        end
        current_shift_idx = current_shift_idx - 1;
        rst = 0;
    end
    if rst < nrst
        ns = length(ev) + nsp;
        rst = rst + 1;
        V = [];
    else
        ns = ns_init;
        rst = 1;
        shift_idx = shift_idx + 1;
        ev = [];
        V = [];
    end
end