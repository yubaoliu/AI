% Program: g_least_pth.m
% Description: Evaluates the gradient of the
% objective function Psi_k in Eq. (8.7). This 
% MATLAB function is required by MATLAB functions 
% least_pth.m and least_pth_m.m.
% ===============================================
function z = g_least_pth(x,pr)
N = pr(1);
p = pr(2);
K = pr(3);
M0 = pr(4:(K+3));
omi = pr((K+4):(2*K+3));
a = x(1:(N+1));
b = x((N+2):end);
bb = [1; b];
M = abs(freqz(a,bb,omi));
err = M - M0;
se = sign(err);
err_abs = abs(err);
em = max(err_abs);
ee = err_abs/em;
cp0 = norm(ee,p)^(1-p);
cp = ee.^(p-1);
z1 = zeros((N+1),1);
z2 = zeros(N,1);
f0 = 0:1:N;
f0 = f0(:);
for i = 1:K,
    fi = omi(i)*f0;
    c = cos(fi);
    s = sin(fi);
    c1 = c(2:(N+1));
    s1 = s(2:(N+1));
    ac = a'*c;
    as = a'*s;
    bc = 1 + b'*c1;
    bs = b'*s1;
    aw = ac^2 + as^2;
    bw = bc^2 + bs^2;
    csmi = cp(i)*se(i)*M(i);
    z1 = z1 + (csmi/aw)*(ac*c + as*s);
    z2 = z2 + (csmi/bw)*(bc*c1 + bs*s1);
end
z = cp0*[z1;-z2];