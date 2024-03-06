function W = LongVec(Prob,MU,U)
%% W = LongVec(Prob,MU,V) returns the long vectors W
%   generated from MU and U.
%           W = [R(MU)*U; V; C(MU)*U]
%   where R() represents rational functions and C()
%   represents Chebyshev polynomails.
% INPUT:
%	MU - shifts, complex values
%	U  - short vectors of length n

%% --- parameters
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials

for k = 1:length(MU)
    mu = MU(k);
    u = U(:,k)/norm(U(:,k)); %--- eigenvector
    w = u; %--- long vector [v1,v2,...,u,u1,u2,...]
    %--- rational part (v1,v2,...)
    for i = np:-1:1
        w = [u/(zk(i)-mu);w];
    end
    %--- Chebyshev part (u1,u2,...)
    if deg > 1
        u0 = u;
        u1 = (mu-c)/r*u;
        w = [w; u1];
    end
    for j = 2:deg-1
        u2 = 2*(mu-c)/r*u1 - u0;
        u0 = u1;
        u1 = u2;
        w = [w; u2];
    end
    %--- store w in W
    W(:,k) = w;
end