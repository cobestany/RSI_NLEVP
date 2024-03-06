%% Test
%% Problem: butterfly problem
%% Method:  Reduced Subspace Iteration

%% Set up Problem
clear
% addpath('data')
%------------------- coefficients matrix
data = load('data/butterfly.mat');
B0 = data.A0;
A0 = data.A1;
A1 = data.A2;
A2 = data.A3;
A3 = data.A4;
%------------------- nonlinear functions
ff1 = @(t) (t.^2);
ff2 = @(t) (t.^3);
ff3 = @(t) (t.^4);
%------------------- original problem (nonlinear part)
F = @(z) (ff1(z)*A1+ff2(z)*A2+ff3(z)*A3);
T = @(z) (B0+A0*z+F(z));
%------------------- relative residual
nrmM = [norm(B0,1) norm(A0,1) norm(A1,1) norm(A2,1) norm(A3,1)];
nrmf = @(z) abs([1 z ff1(z) ff2(z) ff3(z)]);
nrm = @(z) sum(nrmM.*nrmf(z));

%% Parameters
%-------------------- contour data
c = 0;                  %------- center of contour
rin = .5;                 %------- radius of outer circle
r = 1.25*rin;           %------- circle of approx region
np = 0;                %------- #. quadrature points
deg = 4;                 %------- deg of Chebyshev polynomials
ns_init = 30;                %------- subspace dimension
nsp = 15;                %------- if num. ev > ns - nsp, split
ninv = 5;              %------- #. inverse iteration
nref = 0;               %------- #. iterative refinement in shift&invert
nrst = 4;               %------- #. restart in each cluster
m = 1;                  %------- num of clusters
tol = 1e-14;            %------- error tolerance for locking
ns = ns_init;
%-------------------- help variables
n = size(B0,2);         %------- size of the original nlevp
mu0 = c;         %------- common-used shift center

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
%--- init LU for the matrix used in iterative refinement
[LUPQc.L,LUPQc.U,LUPQc.P,LUPQc.Q] = lu(Ci{end});

%% Problem Related data
Prob = struct('Bi',{Bi},'Ci',{Ci},'c',c,'r',r,'sig',zk,'np',np,'deg',deg);

%% Initial Candidate Shifts
[z,~] = modARND(Prob,mu0,40);
[~,shift]= Kmeans(m,z);  %---------- initial shifts
shift_idx = 1;
current_shift_idx = 0;
plotidx = 1;
current_plotidx = 1;

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

%% Reduced Subspace Iteration
eval = [];    %---------- locked eigenvalues
evec = [];    %---------- locked eigenvectors
ev = [];   %---------- new find eigenvalues
V = [];    %---------- basis vectors
Vs = [];   %---------- Schur vectors of locked eigenvectors
rst = 1;   %---------- counter for restarts
while shift_idx <= length(shift)
fprintf('cluster: %d\n',shift_idx)
%--- generate basis V
l = length(ev);
for k = 1:ns
    %---------- init LU factorization of the most lower-right block
    if current_shift_idx ~= shift_idx
        current_shift_idx = shift_idx;
        mu0 = shift(shift_idx);
        [LUPQ0] = initLU(Prob,mu0);
    end
    %---------- use eigenpair (or random vector) for inverse iteration
    if k <= l
        v = U(:,k);
        mu = ev(k);
    else
        v = randn(n,1);
        mu = mu0;
    end
   
    %------------ use current LU or re-compute
    if abs(mu-mu0)<0.01*r
        f = 0; % use pre-computed LU
    else
        f = 1; % compute LU for eigen-pair
    end
    %---------- ninv steps of inverse iteration
    v = ShiInv(Prob,mu,v,ninv,nref,f,LUPQ0,LUPQc);  
    V = mgs(V,k,v);
    fprintf('(%2d/%2d)\n',k,ns);
end
%--- project onto V to reduce problem size
Bip = {};
for k = 1:np
    Bip{k} = V'*(Bi{k}*V);
end
Cip = {sparse(zeros(ns))};
for k = 1:deg+1
    Cip{k} = V'*(Ci{k}*V);
end

%--- solve eigenvalue problem
SubProb = struct('Bi',{Bip},'Ci',{Cip},'c',c,'r',r,'sig',zk,'np',np,'deg',deg);
[ev,U] = solveNLEVP(SubProb);
[ev,idx] = in_cluster(ev,shift,shift_idx);
U = V*U(:,idx);

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

%% Multi-plot Figure (Splitting Cluster Procedure)
    if shift_idx == 1 && rst == 1 && length(shift) == m
        close all
        multifig = figure;
    end
    if plotidx == current_plotidx
    current_plotidx = current_plotidx + 1;
    figure(multifig)
    subplot(2,2,plotidx);
    th = linspace(0,2*pi,100);
    xi = 1.25*r*cos(th)+real(c);
    yi = 1.25*r*sin(th)+imag(c);
    plot(xi,yi,'k--');
    hold on
    eval_true = data.eval;
    plot(real(eval_true),imag(eval_true),'rs','LineWidth',1,'MarkerSize',5)
    plt111 = plot(real(ev),imag(ev),'b+','LineWidth',1,'MarkerSize',8);
    plt222 = plot(real(shift),imag(shift),'*','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','#77AC30');
    plot(real(zk),imag(zk),'ko','LineWidth',1,'MarkerSize',5)
    %--- plot setting
    title(['num. subregion = ' num2str(length(shift))],'FontSize',12)
    axis equal
    xlim([-1 1])
    ylim([-1 1]) 
    end

%% Splitting Clusters
if length(ev) > ns-nsp %--- too many ev in cluster
    fprintf('Too many eigenvalues, add new shifts...\n')
    [~,new_shift]= Kmeans(2,ev); %--- split current cluster into 2
    shift(shift_idx) = new_shift(1); %--- set 1 to be current
    shift = [shift; new_shift(2)]; %--- add 2 to shift set
    current_shift_idx = current_shift_idx - 1;
    rst = 0;
    plotidx = plotidx + 1;
end

if rst < nrst %---- not reaching max num. of restart
    rst = rst + 1;
    V = [];
else %---- reaching nrst, change shift
    rst = 1;
    shift_idx = shift_idx + 1;
    ev = [];
    V = [];
end

end

%% Result
fig3 = figure;
%--- focus region
plot(real(zk),imag(zk),'k.','LineWidth',1,'MarkerSize',18)
hold on
th = linspace(0,2*pi,100);
xi = 1.25*r*cos(th)+real(c);
yi = 1.25*r*sin(th)+imag(c);
% plot(xi,yi,'k--');
xi = rin*cos(th)+real(c);
yi = rin*sin(th)+imag(c);
plot(xi,yi,'r-');
eval_true = data.eval;
plot(real(eval_true),imag(eval_true),'rs','LineWidth',1,'MarkerSize',5)
plot(real(eval),imag(eval),'b+','LineWidth',1,'MarkerSize',8);
plot(real(shift),imag(shift),'*','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','#77AC30');
%--- plot setting
% legend('rational quadrature pts','rational approx region','region of interest','ground truth','RSI Method','shifts','FontSize',14)
legend('region of interest','ground truth','RSI Method','shifts','FontSize',14)
axis equal
xlim([-1 1])
ylim([-1 1])

%--- plot residual of each eigenpair
fig4 = figure;
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
%--- plot residual
semilogy(1:length(res_nrm),res_nrm,'+-');
%--- plot setting
xlabel('index','FontSize',14)
ylabel('Relative Residual','FontSize',14)
grid on
if length(res_nrm)>1, xlim([1 length(res_nrm)]), end
ylim([1e-18 1])