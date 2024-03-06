%% Problem: sandwich beam problem
%% Method:  Reduced Subspace Iteration

%% Set up Problem
close all
addpath('data')
%------------------- coefficients matrix
data = load('sandwich.mat');
Ke = data.Ke;
M = data.M;
Kv = data.Kv;
%------------------- nonlinear functions
ff = data.fun;
%------------------- original problem (nonlinear part)
T = @(z) (Ke-z.^2*M+ff(z)*Kv);
T1 = @(z) (Ke-z.^2*M);
%------------------- relative residual
nrmM = [norm(Ke,1) norm(M,1) norm(Kv,1)];
nrmf = @(z) abs([1 z.^2 ff(z)]);
nrm = @(z) sum(nrmM.*nrmf(z));

%% Parameters
params.deg       = 2;      %--- degree of Chebyshev polynomials
params.n_quad    = 16;     %--- number of quadrature points
params.center    = 5i*1e4; %--- center of region of interest
params.radius    = 1e5;    %--- radius of region of interest
params.tol       = 1e-10;  %--- tolerance of relative residual
params.nit_arnd  = 60;     %--- number of iters for surrogate Arnoldi
params.n_shift   = 1;      %--- number of initial shifts
params.n_restart = 3;      %--- number of restart per subregion
params.nit_inv   = 10;     %--- number of iterations for shift-and-invert
params.dim       = 30;     %--- initial dimension of projection subspace
params.buffer    = 0;      %--- condition for splitting subregions

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
legend('rational quadrature pts','rational approx region','region of interest',...
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
[~,idx] = sort(abs(eval));
res_nrm = res_nrm(idx);
semilogy(1:length(res_nrm),res_nrm,'+-');
xlabel('index','FontSize',14)
ylabel('Relative Residual','FontSize',14)
grid on
if length(res_nrm)>1, xlim([1 length(res_nrm)]), end
ylim([1e-18 1])