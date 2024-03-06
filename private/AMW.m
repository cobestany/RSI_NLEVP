function [AW,MW] = AMW(Prob,W)
%% [AW,MW] = AMW(Prob,W) returns A*W and M*W where A and M form 
%   a generalized eigenvalue problem of size (np+deg)*n:
%               A*w = lam*M*w
%   that approximates the original nonlinear eigenvalue problem.
%   W are given vectors of length (np+deg)*n.


%% --- parameters
Bi = Prob.Bi;  %--- coefficient matrices of rational approx
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
zk = Prob.sig;  %--- quadrature points of rational approx
np = Prob.np;  %--- num. of rational functions
deg = Prob.deg;  %--- degree of Chebyshev polynomials
n = size(Ci{1},1);  %--- problem size
p = size(W,2);   %--- num. of input vectors

%---------- index
%---    rational vectors v = w(1:i1)
%---    base vector u0 = w(i1+1:i2)
%---    Cheb vectors [u1;...;u{deg-2}] = w(i2+1:i3)
%---    last vector u{deg-1} = w(i3+1:end)
i1 = n*np;             %--- end of rational part
i2 = n*np+n;           %--- end of u0
i3 = n*np+n*(deg-1);   %--- end of u{deg-2}

%% --- Part 1 : compute M*w
MW = W;
if deg == 0
    MW(i1+1:end,:) = zeros(n,p);
elseif deg == 1
    MW(i1+1:end,:) = -Ci{deg+1}*MW(i1+1:end,:)/r;
elseif deg == 2
    MW(i1+1:end,:) = MW(i1+1:end,:)/r;
    MW(i2+1:end,:) = -2*(Ci{deg+1}*MW(i2+1:end,:));
else
    MW(i1+1:end,:) = MW(i1+1:end,:)/r;
    MW(i2+1:end,:) = 2*MW(i2+1:end,:);
    MW(i3+1:end,:) = -Ci{deg+1}*MW(i3+1:end,:);
end

%% --- Part 2 : compute A*w
%%----------- rational part
AW = zeros(size(W));
U = W(i1+1:i2,:);
for j = 1:np
    j1 = n*j-n+1;   %--- start of vj
    j2 = n*j;       %--- end of vj
    AW(j1:j2,:) = zk(j)*W(j1:j2,:) + U;
end

%%----------- Chebyshev part
if deg == 0
    AW(i1+1:end,:) = Ci{1}*U;
    for j = 1:np
        j1 = n*j-n+1;   %--- start of vj
        j2 = n*j;       %--- end of vj
        AW(i1+1:end,:) = AW(i1+1:end,:) + Bi{j}*W(j1:j2,:);
    end
elseif deg == 1
    AW(i1+1:end,:) = Ci{1}*U - Ci{2}*U*(c/r);
    for j = 1:np
        j1 = n*j-n+1;   %--- start of vj
        j2 = n*j;       %--- end of vj
        AW(i1+1:end,:) = AW(i1+1:end,:) + Bi{j}*W(j1:j2,:);
    end
elseif deg == 2
    AW(i1+1:i2,:) = c/r*W(i1+1:i2,:) + W(i2+1:end,:);
    AW(i2+1:end,:) = Ci{1}*W(i1+1:i2,:) + Ci{2}*W(i2+1:end,:) ...
        - Ci{3}*(W(i1+1:i2,:) + 2*c/r*W(i2+1:end,:));
    for j = 1:np
        j1 = n*j-n+1;   %--- start of vj
        j2 = n*j;       %--- end of vj
        AW(i2+1:end,:) = AW(i2+1:end,:) + Bi{j}*W(j1:j2,:);
    end
else
    AW(i1+1:i2,:) = c/r*U + W(i2+1:i2+n,:); %--- update u0
    AW(i3+1:end,:) = Ci{1}*U; %--- update u{deg-1}
    for j = 1:deg-2
        j1 = n*(np+j)+1;  %--- start of uj
        j2 = j1+n-1;      %--- end of uj
        U0 = W(j1-n:j1-1,:);
        U1 = W(j1:j2,:);
        U2 = W(j2+1:j2+n,:);
        AW(j1:j2,:) = U0 + 2*c/r*U1 + U2; %--- update uj
        AW(i3+1:end,:) = AW(i3+1:end,:) + Ci{j+1}*U1; %--- update u{deg-1}
    end
    AW(i3+1:end,:) = AW(i3+1:end,:) + Ci{deg}*U2;
    AW(i3+1:end,:) = AW(i3+1:end,:) - Ci{deg+1}*(U1 + 2*c/r*U2);
    for j = 1:np
        j1 = n*j-n+1;   %--- start of v_j
        j2 = n*j;       %--- end of v_j
        AW(i3+1:end,:) = AW(i3+1:end,:) + Bi{j}*W(j1:j2,:);
    end
end