% Program: get_ksi_lambda.m
% Description: Updates parameters lambda and ksi
% by using the foumulas described in Step 4 of 
% Algorithms 8.2 and 8.4. This MATLAB function is 
% required by MATLAB functions charalambous.m and 
% charalambous_m.m.
% ===============================================
function [ksi,lamd] = get_ksi_lambda(x,M0,N,K,ksi,lamd,omi)
a = x(1:(N+1));
b = [1; x((N+2):end)];
M = abs(freqz(a,b,omi));
err = abs(M-M0);
phi = err - ksi;
ind1 = find(phi>0&lamd>0);
i1 = length(ind1);
ind2 = find(phi>0&lamd==0);
i2 = length(ind2);
bphi = 0;
if i1 ~= 0,
    p1 = phi(ind1);
    L1 = lamd(ind1);
    pl1 = p1.*L1;
    bphi = bphi + sum(pl1);
end
if i2 ~= 0,
    p2 = phi(ind2);
    bphi = bphi + sum(p2);
end
L_new = zeros(K,1);
if i1 ~= 0,
   L_new(ind1) = pl1/bphi;
end
if i2 ~= 0,
   L_new(ind2) = p2/bphi;
end
lamd = L_new;
ksi = sum(lamd.*err);