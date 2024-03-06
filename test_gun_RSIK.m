%% Problem: gun problem
%% Method:  Krylov subspace variant
%% NOTE: No truncation Krylov subspace

clear
addpath('data')

%% Set up Problem
clear
addpath('data')
%------------------- coefficients matrix
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
%------------------- original problem (nonlinear part)
F = @(z) (ff1(z)*A1+ff2(z)*A2);
T = @(z) (-B0+A0*z+F(z));
%------------------- relative residual scaling
nrmM = [norm(B0,1) norm(A0,1) norm(A1,1) norm(A2,1)];
nrmf = @(z) abs([1 z ff1(z) ff2(z)]);
nrm = @(z) sum(nrmM.*nrmf(z));

%% Parameters
%-------------------- contour data
c = 1.4e5;              %------- center of contour
rin = 3e4;              %------- radius of contour
r = 1.25*rin;           %------- circle of approx region
np = 32;                %------- #. quadrature points
deg = 10;                 %------- deg of Chebyshev polynomials
ns_init = 30;           %------- dim of subspace
nsp = 10;                %------- if num. ev > ns - nsp, split
nrst = 1;               %------- #. restart in each cluster
m = 5;                  %------- num of clusters
tol = 1e-10;            %------- error tolerance for locking
ns = ns_init;
%-------------------- help variables
n = size(B0,1);         %------- size of the original nlevp
mu0 = c+r*.1i;         %------- common-used shift center

%% Chebyshev approximation
%---------- get Chebyshev coefficient matrices
[Ci] = cheb_approx(T,deg,c,r);

%% Rational approximation
%--- quadrature points
[zk, omega] = contQuad(np,c,1.25*r,0);   %--- quadrature points
%--- rational approx
Bi = {};
for k = 1:np
    z = zk(k);
    %--- original function
    Tz = T(z);
    %--- Chebyshev approximation
    Tcheb = cheb_val(Ci,c,r,z);
    %---------- rational approximation of remaining
    Bi{k} = omega(k)*(Tz-Tcheb);
end

%% Problem Related data
Prob = struct('Bi',{Bi},'Ci',{Ci},'c',c,'r',r,'sig',zk,'np',np,'deg',deg);

%% Initial Candidate Shifts
[z,~] = modARND(Prob,mu0,30);
[~,shift]= Kmeans(m,z);  %---------- initial shifts
shift_idx = 1;
current_shift_idx = 0;

%% Merging Clusters
merge_cluster = 1;
while merge_cluster
    [d,i,j] = minDist(shift);
    if d < 0.1*r
        fprintf('merge clusters...\n')
        new_shift = (shift(i)+shift(j))/2;
        shift([i j]) = [];
        shift = [shift; new_shift];
        m = m-1;
    else
        merge_cluster = 0;
    end     
end


%% Krylov Subspace Variant (Full Length Basis)
eval = [];    %---------- locked eigenvalues
evec = [];    %---------- locked eigenvectors
ev = [];   %---------- new eigenvalues
V = [];    %---------- basis vectors
Vs = [];   %---------- Schur vectors of locked eigenvectors
rst = 1;   %---------- counter for restarts
while shift_idx <= length(shift)
fprintf('cluster: %d\n',shift_idx)
%--- generate basis V
l = length(ev);
W = randn((np+deg)*n,1);
W = W/norm(W);
for k = 2:ns
    %---------- LU factorization of the Schur complement
    if current_shift_idx ~= shift_idx
        current_shift_idx = shift_idx;
        mu0 = shift(shift_idx);
        [LUPQ0] = initLU(Prob,mu0);
        fprintf('( 1/%2d)\n',ns)
    end
    %------------ generate next Krylov basis vector
    w = W(:,k-1);
    [w,~] = ShiInv_solve(Prob,mu0,w,LUPQ0,0);
    %---------- truncate Krylov subspace (w to v)
    W = mgs(W,k,w);
    fprintf('(%2d/%2d)\n',k,ns);
end

%--- project onto W to reduce problem size
[AW,MW] = AMW(Prob,W);
Ap = W'*AW;
Mp = W'*MW;

%--- solve eigenvalue problem
[U, ev] = eig(full(Ap),full(Mp));
ev = diag(ev);

%---------- exclude outsiders
idx = abs(ev-c)<rin;
U = W(np*n+1:np*n+n,:)*U(:,idx);
ev = ev(idx);

%---------- check incluster
[ev,idx] = in_cluster(ev,shift,shift_idx);
U = U(:,idx);

%% Residual
%--- relative residual
% [res_vec,res_nrm] = get_res(Prob,ev,U);
% for i = 1:length(ev)
%     u = U(:,i);
%     res_nrm(i) = norm(res_vec)/norm(u)/nrm(ev(i));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res_nrm = [];
for i = 1:length(ev)
    u = U(:,i);
    res_vec = T(ev(i))*u/norm(u);
    res_nrm(i) = norm(res_vec)/nrm(ev(i));
end

%--- locking
idx = find(res_nrm < tol);
for i = 1:length(idx)
    ui = U(:,idx(i));
    d = 1.0; %--- distance to locked Schur vectors
    if ~isempty(Vs)
        d = norm(ui - Vs*(Vs'*ui))/norm(ui);
    end
    if d < 1e-3,continue,end %--- angle not large enough
    Vs = mgs(Vs,size(Vs,2)+1,ui);
    eval = [eval; ev(idx(i))];
    evec = [evec, ui];
end

%% Splitting Clusters
if length(ev) > ns_init-nsp %--- too many ev in one cluster
    [~,new_shift]= Kmeans(2,ev); %--- split current cluster into 2
    shift(shift_idx) = new_shift(1); %--- set 1 to be current
    shift = [shift; new_shift(2)]; %--- add 2 to shift set
    current_shift_idx = current_shift_idx - 1;
    rst = 0;
end

if rst < nrst %---- not reaching max num. of restart
    ns = length(ev) + nsp;
    rst = rst + 1;
    V = [];
else %---- reaching nrst, change shift
    ns = ns_init;
    rst = 1;
    shift_idx = shift_idx + 1;
    ev = [];
    V = [];
end

pause(1)
end

%% Result
close all
figure(1)
%--- focus region
plot(real(zk),imag(zk),'k.','LineWidth',1,'MarkerSize',18)
hold on
th = linspace(0,2*pi,100);
xi = 1.25*r*cos(th)+real(c);
yi = 1.25*r*sin(th)+imag(c);
plot(xi,yi,'k--');
xi = rin*cos(th)+real(c);
yi = rin*sin(th)+imag(c);
plot(xi,yi,'r-');
plot(real(z),imag(z),'md')
plot(real(eval),imag(eval),'b+','LineWidth',1,'MarkerSize',8);
plot(real(shift),imag(shift),'g*','LineWidth',1,'MarkerSize',8);
%--- plot setting
legend('rational quadrature pts','rational approx region','region of interest',...
    'Modified Arnoldi Method','RSIK Method','shifts','FontSize',14)
axis equal

%--- plot residual of each eigenpair
figure(2)
%--- relative residual
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
%--- plot setting
xlabel('index','FontSize',14)
ylabel('Relative Residual','FontSize',14)
grid on
if length(res_nrm)>1, xlim([1 length(res_nrm)]), end
ylim([1e-18 1])