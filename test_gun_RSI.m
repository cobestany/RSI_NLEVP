%% Problem: gun problem
%% Method:  Reduced Subspace Iteration

%% Set up Problem
close all
addpath('data')
%------------------- coefficient matrices
data = load('gun.mat');
B0 = -data.K;
A0 = -data.M;
A1 = data.W1;
A2 = data.W2;
%------------------- nonlinear functions
sig1 = 0;
sig2 = 108.8774;
ff1 = @(t) (1i*sqrt(t - sig1^2));
ff2 = @(t) (1i*sqrt(t - sig2^2));
%------------------- original problem
F = @(z) (ff1(z)*A1+ff2(z)*A2);
T = @(z) (-B0+A0*z+F(z));
%------------------- relative residual scaling function
nrmM = [norm(B0,1) norm(A0,1) norm(A1,1) norm(A2,1)];
nrmf = @(z) abs([1 z ff1(z) ff2(z)]);
nrm = @(z) sum(nrmM.*nrmf(z));

%% Parameters
params.deg       = 5;    %--- degree of Chebyshev polynomials
params.n_quad    = 32;    %--- number of quadrature points
params.center    = 1.4e5; %--- center of region of interest
params.radius    = 3e4;   %--- radius of region of interest
params.tol       = 1e-10; %--- tolerance of relative residual
params.nit_arnd  = 30;    %--- number of iters for surrogate Arnoldi
params.n_shift   = 5;     %--- number of initial shifts
params.n_restart = 2;     %--- number of restart per subregion
params.nit_inv   = 5;     %--- number of iterations for shift-and-invert
params.dim       = 10;    %--- initial dimension of projection subspace
params.buffer    = 2;     %--- condition for splitting subregions

%% Reduced Subspace Iteration
[eval,evec,Prob,shift] = RSI(T,params,nrm);

%% Result
zk = Prob.sig;
r = Prob.r;
c = Prob.c;

figure(1)
plot(real(zk),imag(zk),'k.','LineWidth',1,'MarkerSize',18)
hold on
th = linspace(0,2*pi,100);
xi = 1.25*r*cos(th)+real(c);
yi = 1.25*r*sin(th)+imag(c);
plot(xi,yi,'k--');
xi = params.radius*cos(th)+real(c);
yi = params.radius*sin(th)+imag(c);
plot(xi,yi,'r-');
plot(real(eval),imag(eval),'b+','LineWidth',1,'MarkerSize',8);
plot(real(shift),imag(shift),'*','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','#77AC30');
legend('rational quadrature pts','region of approx','region of interest',...
    'RSI Method','shifts','FontSize',14)
axis equal

%--- relative residual
figure(2)
res_nrm = [];
for i = 1:length(eval)
    u = evec(:,i);
    res_vec = T(eval(i))*u/norm(u);
    res_nrm(i) = norm(res_vec)/nrm(eval(i));
end
%--- sort res_nrm w.r.t. distance from c
[~,idx] = sort(abs(eval-c));
res_nrm = res_nrm(idx);
semilogy(1:length(res_nrm),res_nrm,'+-');
xlabel('index','FontSize',14)
ylabel('Relative Residual','FontSize',14)
grid on
if length(res_nrm)>1, xlim([1 length(res_nrm)]), end
ylim([1e-18 1])