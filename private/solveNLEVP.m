function [eval,U] = solveNLEVP(Prob)
%% [eval,U] = solveNLEVP(Prob) returns the eigenvalues
%   and eigenvectors of a rational-Chebyshev approximate
%   problem:
%               T(z) ~ R(z) + C(z)
%   Linearizing the above problem into a large generalized
%   eigenvalue problem:
%               A*w = lam*M*w
%   which has the same spectrum as the NLEVP.

%% --- parameters
Ci = Prob.Ci;  %--- coefficient matrices of Chebyshev approx
c = Prob.c;    %--- center of the contour
r = Prob.r;    %--- radius of the contour
np = Prob.np;  %--- num. of rational functions
n = size(Ci{1},2);

%--- generate A and M
[A,M] = genAM(Prob);

%--- solve generalized eigenvalue problem
[U, eval] = eig(full(A),full(M));
eval = diag(eval);

%--- exclude points outside of region of interest
idx = abs(eval-c)<0.8*r;
U = U(np*n+1:np*n+n,idx);
eval = eval(idx);