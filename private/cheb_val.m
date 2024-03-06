   function [Tcheb] = cheb_val(Ci, c, r, z)
%% [Tcheb] = cheb_val(Ci, c, r, z) returns the matrix
%   value of C(z) such that
%             C(z) = C0*t0(z) + C1*t1(z) + ...
%   where tj(z) is the Chebyshev polynomial of degree
%   j of the first kind.
% INPUT:
%    Ci   = cell of Chebyshev coefficient matrices
%    c    = center of contour
%    r    = radius of contour
%    z    = approximate position

%--- degree
deg = length(Ci)-1;

%--- transform z onto a unit circle
zc = (z-c)/r;

%--- case: degree 0
if deg == 0
    Tcheb = Ci{1};
    return
end

%--- Chebyshev approximation
%---------- degree 1 part
t0 = 1;
t1 = zc;
Tcheb = Ci{1} + Ci{2}*t1;
%---------- degree 2:deg
for ki = 2:deg
    t2 = 2*zc*t1 - t0;
    t0 = t1;
    t1 = t2;
    Tcheb = Tcheb + Ci{ki+1}*t2;
end