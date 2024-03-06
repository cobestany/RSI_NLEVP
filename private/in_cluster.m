function [y,f] = in_cluster(x,shift,idx)
%%-------- x is a set of points
%%-------- shift is a set of shifts
%%-------- idx is an integer indicator of current shift
%%-------- f indicates whether x belongs to shift(idx)
%%-------- y is the of points in current cluster

m = length(x);

for i = 1:m
    current_x = x(i);
    [~,id] = min(abs(current_x-shift));
    if id == idx
        f(i) = true;
    else
        f(i) = false;
    end
end

y = x(f);