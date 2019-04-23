% Program: bfgs_charalambous.m
% Description: Implements Step 3 of Algorithm 8.2 and
% Algorithm 8.4 by using the quasi-Newton algorithm 
% with BFGS updating formula. This MATKAB function is 
% required for the implementation of MATLAB functions 
% charalambous.m and charalambous_m.m.
% ====================================================
function x = bfgs_charalambous(x_ini,N,K,ksi,M0,lamd,omi,epsi2)
xk = x_ini(:);
pr = [N;K;ksi;M0;lamd;omi];
I = eye(length(xk));
Sk = I;
gk = g_charalambous(xk,pr);
fk = f_charalambous(xk,pr);
dk = -Sk*gk;
ak = inex_lsearch(xk,dk,'f_charalambous','g_charalambous',pr);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = f_charalambous(xk_new,pr);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
while err >= epsi2,
      gk_new = g_charalambous(xk_new,pr);
      gmk = gk_new - gk;
      D = dtk'*gmk;
      if D <= 0,
         Sk = I;
      else
         sg = Sk*gmk;
         sw0 = (1+(gmk'*sg)/D)/D;
         sw1 = dtk*dtk';
         sw2 = sg*dtk';
         Sk = Sk + sw0*sw1 - (sw2'+sw2)/D;
      end
      fk = fk_new;
      gk = gk_new;
      xk = xk_new;
      dk = -Sk*gk;
      ak = inex_lsearch(xk,dk,'f_charalambous','g_charalambous',pr);
      dtk = ak*dk;
      xk_new = xk + dtk;
      fk_new = f_charalambous(xk_new,pr);
      dfk = abs(fk - fk_new);
      err = max(dfk,norm(dtk));
end
x = xk_new;