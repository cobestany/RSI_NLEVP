% a is an array, minDist(a) returns the minimum
% distance between any two different entries in a
function [d,i,j] = minDist(a)
n = length(a);
[X,Y] = meshgrid(a);
Z = triu(abs(X-Y),1) + tril(inf*ones(n));
d = min(min(Z));
[i,j] = find(Z==d);