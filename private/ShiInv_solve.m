function [w,LUPQ] = ShiInv_solve(Prob,mu,w,LUPQ,f)
%% [w,LL,UU] = ShiInv_solve(Prob,mu,w,mu,LL,UU,f) returns
%               w = (A - mu*M)\(M*w)
%   where A and M are formed implicitly.


%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials
n = size(Ci{1},1);  %--- problem size

LL = LUPQ.L;
UU = LUPQ.U;
PP = LUPQ.P;
QQ = LUPQ.Q;

%---------- index
%---    rational vectors v = w(1:i1)
%---    base vector u0 = w(i1+1:i2)
%---    Cheb vectors [u1;...;u{deg-2}] = w(i2+1:i3)
%---    last vector u{deg-1} = w(i3+1:end)
i1 = n*np;            %--- end of rational part
i2 = n*np+n;          %--- end of u0
i3 = n*np+n*max(deg-1,0);   %--- end of u{deg-2}

%% --- Step 1 : w <- M*w
if deg == 0
    w(i1+1:end) = zeros(n,1);
elseif deg == 1
    w(i1+1:end) = -Ci{deg+1}*w(i1+1:end)/r;
elseif deg == 2
    w(i1+1:end) = w(i1+1:end)/r;
    w(i2+1:end) = -2*(Ci{deg+1}*w(i2+1:end));
else
    w(i1+1:end) = w(i1+1:end)/r;
    w(i2+1:end) = 2*w(i2+1:end);
    w(i3+1:end) = -Ci{deg+1}*w(i3+1:end);
end

%% --- Step 2 : w <- L\w
for j = 1:np
    j1 = n*j-n+1;   %--- start of v_j
    j2 = n*j;       %--- end of v_j
    w(i3+1:end) = w(i3+1:end) - Bi{j}*w(j1:j2)/(zk(j)-mu);
end

%% --- Step 3 : w <- U\w
%% ---------- construct S_mu
%---------- replace C_0 <- C_0 - B/(D-mu*I)*E
C0 = Ci{1};
for j = 1:np
    C0 = C0 - Bi{j}/(zk(j)-mu);
end
%% ---------- solve S_mu*x = t
if deg == 0
    if f
        [LL,UU,PP,QQ] = lu(C0);
    end
    w(i1+1:i2) = QQ*(UU\(LL\(PP*w(i1+1:i2))));
elseif deg == 1
    if f
        alpha = (c-mu)/r;
        [LL,UU,PP,QQ] = lu(C0-alpha*Ci{2});
    end
    w(i1+1:i2) = QQ*(UU\(LL\(PP*w(i1+1:i2))));
else
    alpha = 2*(c-mu)/r; %--- original diagonal entries
    beta = 2/alpha;     %--- beta_0
    %---------- u1
    w(i1+1:i2) = beta(1)*w(i1+1:i2);
    %---------- Gaussian elimination
    for j = 2:deg-1
        %---------- index
        j1 = i1 + (j-2)*n;
        j2 = j1 + n;
        %---------- update uj
        w(j2+1:j2+n) = w(j2+1:j2+n) - w(j1+1:j1+n);
        beta(j) = 1/(alpha - beta(j-1));
        w(j2+1:j2+n) = beta(j)*w(j2+1:j2+n);
        %---------- update u{end}
        w(i3+1:i3+n) = w(i3+1:i3+n) - C0*w(j1+1:j1+n);
        C1 = Ci{j} - beta(j-1)*C0;
        C0 = C1;
    end
    %---------- update C_{end-1}
    C1 = C0 - Ci{deg+1};
    C0 = C1;
    %---------- update C_{end}
    w(i3+1:end) = w(i3+1:end) - C0*w(i3-n+1:i3);
    C1 = Ci{deg} - alpha*Ci{deg+1} - beta(deg-1)*C0;
    %---------- solve x_{nk-1}
    if f
        [LL,UU,PP,QQ] = lu(C1);
    end
    w(i3+1:end) = QQ*(UU\(LL\(PP*w(i3+1:end))));
    %---------- backward solve
    for j = deg-1:-1:1
        j1 = i1 + (j-1)*n;
        j2 = j1 + n;
        w(j1+1:j1+n) = w(j1+1:j1+n) - beta(j)*w(j2+1:j2+n);
    end
end

%% ---------- v <- (D-mu*I)\(v-E*t)
ue = w(i1+1:i2);
for j = 1:np
    j1 = n*j-n+1;   %--- start of v_j
    j2 = n*j;       %--- end of v_j
    w(j1:j2) = (w(j1:j2)-ue)/(zk(j)-mu);
end

LUPQ.L = LL;
LUPQ.U = UU;
LUPQ.P = PP;
LUPQ.Q = QQ;