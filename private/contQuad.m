function [zk,omega] = contQuad(np,c,r,method)
% [zk,omega] = contQuad(n,c,r,method) generates quadrature points and 
%   weights on a defined circle. np is the number of points. c and 
%   r are center and radius of the circle. method describes quadrature rule: 
%   0 = Gaussian, 1 = trapezoidal, else = mid-point.

if np == 0
    zk = [];
    omega = [];
    return
end

if( method == 0 )
%-------------------- Gauss-Legendre    
    disp 'Use Gauss Legendre';
    np = fix(np/2); % num of points on half circle
    np = np - 1 ; % n + 1 points without this line
    beta = .5 ./ sqrt(1 - (2 * ( 1 : np ) ).^( - 2 ) );
    T = diag( beta, 1 ) + diag( beta, - 1 ) ;
    [ V, D ] = eig( T ) ;
    x = diag( D ); 
    [ x, I ] = sort( x );
    omega = ( V( 1, I ).^2 )';
    theta = (pi / 2) .* (1 - x);
    %---------------- points on unit half-circle
    zk = exp(1i*theta);
    omega = -omega.*zk;
    %---------------- translate to defined circle
    zk = [zk; conj(zk(end:-1:1))];
    zk = c + r*zk;
    omega = [omega; conj(omega(end:-1:1))];
    omega = 0.5*r*omega;
elseif( method == 1 )
%-------------------- trapezoidal rule
    disp 'Use trapezoidal rule';
    x = linspace(-1,1,np);
    theta = pi.*x';
    omega = 2*ones(np,1)/np;
    omega(1) = 1/np;
    omega(end) = 1/np;
    %---------------- points on unit half-circle
    zk = exp(1i*theta);
    omega = -omega.*zk;
    %---------------- translate to defined circle
    zk = c + r*zk;
    omega = 0.5*r*omega;
else
%-------------------- mid-point
    disp 'Use mid-point';
    x = linspace(-1,1,np);
    theta = pi.*x';
    omega = ones(np,1)/np;
    %---------------- points on unit circle
    zk = exp(1i*theta);
    omega = -omega.*zk;
    %---------------- translate to defined circle
    zk = c + r*zk;
    omega = r*omega;
end
