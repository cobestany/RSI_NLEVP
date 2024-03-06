function [X] =  mgs(X,j,V) 
% [X] = mgs(X,j,v) returns the modified Gran-Schmidt basis.
%   Orthonormalize columns of V against those of X(:,1:j-1) and
%   put in column X(:,j:end). 

   [m,nX] = size(X);
   [m1, nv] = size(V);
   if (m == 0)
     X = zeros(m1,1);
   else
     if (m1 ~= m)
       error('dimensions of X and v do not match');
     end     
   end
%%-------------------- main loop     
   for k = 1:nv
     w = V(:,k);
     for i = 1:j+k-2
       t = X(:,i)'*w;
       w = w - t*X(:,i);
     end
     t = norm(w); 
     X(:,j+k-1) = w/t; 
   end
 end
