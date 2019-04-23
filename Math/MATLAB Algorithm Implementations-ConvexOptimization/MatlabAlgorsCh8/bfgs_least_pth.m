% Program: bfgs_least_pth.m
% Description: Implements Step 3 of Algorithm 8.1 and
% Algorithm 8.3 by using the quasi-Newton algorithm 
% with BFGS updating formula. This MATLAB function is 
% required for the implementation of MATLAB functions 
% least_pth.m and least_pth_m.m.
% ====================================================
function x = bfgs_least_pth(x_ini,N,p,K,M0,omi,epsi2)
xk = x_ini(:);
pr = [N;p;K;M0;omi];
I = eye(length(xk));
Sk = I;
gk = g_least_pth(xk,pr);
fk = f_least_pth(xk,pr);
dk = -Sk*gk;
ak = inex_lsearch(xk,dk,'f_least_pth','g_least_pth',pr);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = f_least_pth(xk_new,pr);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
while err >= epsi2,
      gk_new = g_least_pth(xk_new,pr);
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
      ak = inex_lsearch(xk,dk,'f_least_pth','g_least_pth',pr);
      dtk = ak*dk;
      xk_new = xk + dtk;
      fk_new = f_least_pth(xk_new,pr);
      dfk = abs(fk - fk_new);
      err = max(dfk,norm(dtk));
end
x = xk_new;
