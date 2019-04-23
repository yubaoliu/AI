% Program: f_least_pth.m
% Description: Evaluates the objective function 
% Psi_k in Eq. (8.6). This MATLAB function is 
% required by MATLAB functions least_pth.m and 
% least_pth_m.m.
% =============================================
function z = f_least_pth(x,pr)
N = pr(1);
p = pr(2);
K = pr(3);
M0 = pr(4:(K+3));
omi = pr((K+4):(2*K+3));
a = x(1:(N+1));
b = [1; x((N+2):end)];
M = abs(freqz(a,b,omi));
ae = abs(M - M0);
em = max(ae);
ee = ae/em;
z = em*norm(ee,p);