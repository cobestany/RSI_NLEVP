function [cluster,centroid] = Kmeans(k,P)
%% [cluster,centroid] = Kmeans(k,P) assigns points P to
%   k clusters and returns the centroids of clusters.
% INPUT:
%   k        - number of desired clusters
%   P        - m complex points
% OUTPUT:
%   cluster  - indicator of clusters for each point, m-by-1 array
%   centroid - centroids of each cluster, k-by-1 array
% WARNING:
%   TO FIX: centroid can have NaN

%--------- parameters
m = length(P);
k = min(k,m);
%--------- choose k unique random indices from P
idx = randperm(m,k);
%--------- set initial centroids
centroid = P(idx);

%--------- init cluster array
cluster = zeros(m,1); % cluster indicator of each point
minDist = zeros(m,1); % distance from centre of each point

%--------- init previous cluster array clusterPrev (for stopping criterion)
clusterPrev = cluster;

%--------- init stopping flag
if k > 0
  f = 1;
else
  f = 0; % no input
end

while f
    %--------- for each data point 
    for idxP = 1:m
        %--------- init distance array
        dist = zeros(k,1);
        %--------- compute distance to each centroid
        for idxC=1:k
            dist(idxC) = abs(P(idxP)-centroid(idxC));
        end
        %--------- find index of closest centroid
        [minDistP, clusterP] = min(dist);
        cluster(idxP) = clusterP;
        minDist(idxP) = minDistP;
    end
        
    %--------- for every cluster compute new centroid
    for idxC = 1:k
        %--------- find the points in cluster idxC and compute mean
        centroid(idxC) = mean(P(cluster==idxC));
    end
    
    %--------- Checking for stopping criterion: Clusters do not chnage
    if clusterPrev==cluster
        f = 0;
    else
    %--------- update previous cluster clusterPrev
      clusterPrev = cluster;
    end
end