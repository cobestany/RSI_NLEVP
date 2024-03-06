function y = ShiInv_solve_rev(Prob,mu,x,LUPQ)
%% [w,LL,UU] = ShiInv_solve_rev(Prob,mu,x,LL,UU) returns
%               y = M\(A-mu*M)*x
%   where A and M are formed implicitly. This is the reverse
%   step of shift-and-invert.


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

%% ----- compute x <- (A-mu*M)*x
y = zeros(size(x));
%----- rational part
for j = 1:np
    j1 = n*j-n+1;   %--- start of v_j
    j2 = n*j;       %--- end of v_j
    y(j1:j2) = (zk(j)-mu)*x(j1:j2) + x(i1+1:i1+n);
    y(i3+1:i3+n) = y(i3+1:i3+n) + Bi{j}*x(j1:j2);
end
%----- u part
x0 = x(i1+1:i1+n);
if deg == 0
    y(i1+1:i1+n) = y(i1+1:i1+n) + Ci{1}*x0;
elseif deg == 1
    y(i1+1:i1+n) = y(i1+1:i1+n) + Ci{1}*x0 - Ci{2}*x0*(c-mu)/r;
elseif deg == 2
    y(i1+1:i2) = (c-mu)/r*x0 + x(i2+1:end);
    y(i2+1:end) = y(i2+1:end) + Ci{1}*x0 + Ci{2}*x(i2+1:end) ...
        - Ci{3}*(x0 + 2*(c-mu)/r*x(i2+1:end));
else
    alpha = 2*(c-mu)/r;
    y(i1+1:i1+n) = alpha/2*x0 + x(i2+1:i2+n);
    y(i3+1:i3+n) = y(i3+1:i3+n) + Ci{1}*x0;
    %----- Chebyshev part
    for j = 1:deg-2
        j1 = i1 + (j-1)*n;
        j2 = j1 + n;
        j3 = j2 + n;
        y(j2+1:j2+n) = x(j1+1:j1+n) + alpha*x(j2+1:j2+n) + x(j3+1:j3+n);
        y(i3+1:i3+n) = y(i3+1:i3+n) + Ci{j+1}*x(j2+1:j2+n);
    end
    %----- last row
    y(i3+1:i3+n) = y(i3+1:i3+n) - Ci{deg+1}*x(j2+1:j2+n)...
        + Ci{deg}*x(i3+1:i3+n) - alpha*Ci{deg+1}*x(i3+1:i3+n);
end

%% ----- compute x <- M\x
if deg == 0
    return
elseif deg == 1
    y(i1+1:end) = -QQ*(UU\(LL\(PP*y(i1+1:end))))*r;
elseif deg == 2
    y(i1+1:end) = r*y(i1+1:end);
    y(i2+1:end) = -QQ*(UU\(LL\(PP*y(i2+1:end))))/2;
else
    y(i1+1:end) = r*y(i1+1:end);
    y(i2+1:end) = y(i2+1:end)/2;
    y(i3+1:end) = -QQ*(UU\(LL\(PP*y(i3+1:end))));
end