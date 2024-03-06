function [Ci] = cheb_approx(T,deg,c,r)
%% Ci = cheb_approx(T,deg,c,r,n) returns a cell of coefficient matrices
%   of the Chebyshev approximation on a shifted-and-scaled circle. The
%   Chebyshev approximation of the nonlinear eigenvalue problem will be:
%              T(z) ~ Ci{1}*t0(z) + ... + Ci{deg+1}*t{deg}(z)
%   where tj(z) is the chebyshev polynomial of degree j of the first kind.
%
% INPUT:
%   T   - matrix-valued function of original problem
%   deg - degree of Chebyshev polynomial
%   c   - center of contour
%   r   - radius of contour

%% interpolation
Ci = cell(deg+1,1);
%%%---- num. of interpolation nodes
npts = 2*deg + 10;
%%%---- interpolation nodes
th = [1:2:2*npts-1]*(pi/(2*npts)) ;
xi = cos(th(:));
zi = r*xi + c;
%%%---- matrix-valued function T(z)
Tk = cell(npts,1);
for k = 1:npts
    Tk{k} = sparse(T(zi(k)));
end
[m,n] = size(Tk{1});

%% Compute coefficient matrices
%%%---- 
for k = 0:deg
    zk = cos(k*th(:));
    Tc = sparse(m,n);
    
    for j = 1:npts
        Tc = Tc + zk(j)*Tk{j};
    end
    
    Ci{k+1} = 2.0/npts*Tc;
end
Ci{1} = Ci{1}*0.5;