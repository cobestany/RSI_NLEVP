function [LUPQ] = initLU(Prob,mu)
%% [LL,UU] = initLU(Prob,mu) computes the LU of the Schur
%   complement at shift mu. The purpose is to avoid 
%   repeating LU if mu doesn't change significantly.
%   It follows the same rule as the shift-and-invert: 
%               (A-mu*M)\M*w
%   where only w changes every time.
% INPUT:
%   mu - shift, complex value

%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials

%% --- construct Schur complement
%---------- replace C_0 <- C_0 - B/(D-mu*I)*E
C0 = Ci{1};
for j = 1:np
    C0 = C0 - Bi{j}/(zk(j)-mu);
end
%---------- compute Schur complement S_mu
if deg == 0
    [LL,UU,PP,QQ] = lu(C0);
elseif deg == 1
    alpha = (c-mu)/r;
    [LL,UU,PP,QQ] = lu(C0-alpha*Ci{2});
else
    alpha = 2*(c-mu)/r; %--- original diagonal entries
    beta = 2/alpha;     %--- beta_0
    %---------- Gaussian elimination
    for j = 2:deg-1
        beta(j) = 1/(alpha - beta(j-1));
        C1 = Ci{j} - beta(j-1)*C0;
        C0 = C1;
    end
    %---------- update Ci{end-1}
    C1 = C0 - Ci{deg+1};
    C0 = C1;
    %---------- update Ci{end}
    C1 = Ci{deg} - alpha*Ci{deg+1} - beta(deg-1)*C0;
    %---------- solve x_{nk-1}
    [LL,UU,PP,QQ] = lu(C1);
end

LUPQ.L = LL;
LUPQ.U = UU;
LUPQ.P = PP;
LUPQ.Q = QQ;