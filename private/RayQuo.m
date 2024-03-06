function lam = RayQuo(Prob,u)
%% mu = RayQuo(Prob,u) returns estimated eigenvalue
%   of a given u based on Rayleigh quotient. The idea 
%   is to find the optimal mu such that
%               u'*T(lam)*u/(u'*u) = 0

%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials

%--- construct 1-D sub-problem u'*T(mu)*u/(u'*u) = 0
dn = u'*u;
beta = {};
gamma = {};
for i = 1:np
    beta{i} = (u'*(Bi{i}*u))/dn;
end
for i = 1:deg+1
    gamma{i} = (u'*(Ci{i}*u))/dn;
end
SubProb = struct('Bi',{beta},'Ci',{gamma},'c',c,'r',r,'sig',zk,'np',np,'deg',deg);

%--- linearize the sub-problem
[A,M] = genAM(SubProb);

%--- Rayleigh Quotient
ritz = eig(full(A),full(M));
ritz = ritz(abs(ritz-c)<r);
if isempty(ritz) %--- if optimal mu is out of range
                 %--- return center plus noise
    lam = c+r*1e-5*randn(1,'like',1i);
    return
end
%--- compute residual of each pair (ritz,u)
U = kron(ones(1,length(ritz)),u);
[~,res_nrm] = get_res(Prob,ritz,U);
[~,idx] = min(res_nrm);
%--- select mu with smallest residual
lam = ritz(idx);