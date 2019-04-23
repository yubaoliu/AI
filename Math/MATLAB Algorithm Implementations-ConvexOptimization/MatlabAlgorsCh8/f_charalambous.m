% Program: f_charalambous.m
% Description: Evaluates objective function 
% Psi in Eq. (8.8). This MATLAB function is
% required by MATLAB functions charalambous.m
% and charalambous_m.m.
% ===========================================
function f = f_charalambous(x,pr)
N = pr(1);
K = pr(2);
ksi = pr(3);
M0 = pr(4:(K+3));
lamd = pr((K+4):(2*K+3));
omi = pr((2*K+4):end);
a = x(1:(N+1));
b = [1; x((N+2):end)];
M = abs(freqz(a,b,omi));
phi = abs(M-M0) - ksi;
ind1 = find(phi>0&lamd>0);
i1 = length(ind1);
ind2 = find(phi>0&lamd==0);
i2 = length(ind2);
f = 0;
if i1 ~= 0,
    p1 = phi(ind1).^2;
    L1 = lamd(ind1);
    f = f + 0.5*sum(p1.*L1);
end
if i2 ~= 0,
    f = f + 0.5*sum(phi(ind2).^2);
end