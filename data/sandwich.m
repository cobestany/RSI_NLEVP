%% Problem: sandwich beam problem
%% generate coefficients and functions
function sandwich()
    [coeffs] = nlevp('sandwich_beam');
    Ke = coeffs{1};
    M = coeffs{2};
    Kv = coeffs{3};
    
    alpha = .675;
    tau = 8.230e-9;
    G0   = 3.504e5;
    Ginf = 3.062e9;
    fun = @(lam) (G0 + Ginf*(1i*lam*tau).^alpha)./(1 + (1i*lam*tau).^alpha);
    
    save('sandwich.mat','Ke','M','Kv','fun')
end