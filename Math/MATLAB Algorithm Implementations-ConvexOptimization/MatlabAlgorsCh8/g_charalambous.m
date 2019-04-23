% Program: g_charalambous.m
% Description: Evaluates the gradient of the
% objective function Psi in Eq. (8.11). This 
% MATLAB function is required by MATLAB functions 
% charalambous.m and charalambous_m.m.
% ===============================================
function g = g_charalambous(x,pr)
N = pr(1);
K = pr(2);
ksi = pr(3);
M0 = pr(4:(K+3));
lamd = pr((K+4):(2*K+3));
omi = pr((2*K+4):end);
a = x(1:(N+1));
b = x((N+2):end);
f0 = 0:1:N;
f0 = f0(:);
bb = [1; b];
M = abs(freqz(a,bb,omi));
err = M - M0;
se = sign(err);
err_abs = abs(err);
phi = err_abs - ksi;
ind1 = find(phi>0&lamd>0);
i1 = length(ind1);
ind2 = find(phi>0&lamd==0);
i2 = length(ind2);
g1 = zeros((N+1),1);
g2 = zeros(N,1);
if i1 > 0,
    p1 = phi(ind1);
    L1 = lamd(ind1);
    pl1 = p1.*L1;
    se1 = se(ind1);
    omi1 = omi(ind1);
    M1 = M(ind1);
    for i = 1:i1,
        fi = omi1(i)*f0;
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
        g1 = g1 + (pl1(i)*se1(i)*M1(i)/aw)*(ac*c + as*s);
        g2 = g2 + (pl1(i)*se1(i)*M1(i)/bw)*(bc*c1 + bs*s1);
    end
end
if i2 > 0,
    p2 = phi(ind2);
    se2 = se(ind2);
    omi2 = omi(ind2);
    M2 = M(ind2);
    for i = 1:i2,
        fi = omi2(i)*f0;
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
        g1 = g1 + (p2(i)*se2(i)*M2(i)/aw)*(ac*c + as*s);
        g2 = g2 + (p2(i)*se2(i)*M2(i)/bw)*(bc*c1 + bs*s1);
    end
end
g = [g1;-g2];