% Program: rank_one.m
% Title: Quasi_Newton Rank-One algorithm 
% Description: Implements the quasi-Newton algorithm
% with the rank-one updating formula in Eq. (7.20) as
% described in Algorithm 7.2.
% Theory: See Practical Optimization Sec. 7.4.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   x0: initial point
%   epsi: termination tolerance
% Output:
%   xs: solution point
%   fs: value of objective function at point xs.
%   k: number of iterations required
% Example:
% Find the minimum of the Himmelblau function
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination 
% tolerance epsi = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = rank_one('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ========================================================
function [xs,fs,k] = rank_one(fname,gname,x0,epsi)
disp(' ')
disp('Program rank_one.m')
n = length(x0);
I = eye(n);
k = 1;
xk = x0;
Sk = I;
fk = feval(fname,xk);
gk = feval(gname,xk);
dk = -Sk*gk;
ak = inex_lsearch(xk,dk,fname,gname);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = feval(fname,xk_new);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
while err >= epsi,
      gk_new = feval(gname,xk_new);
      gmk = gk_new - gk;
      ek = dtk - Sk*gmk;
      Sk = Sk + (ek*ek')/(gmk'*ek);
      fk = fk_new;
      gk = gk_new;
      xk = xk_new;
      dk = -Sk*gk;
      ak = inex_lsearch(xk,dk,fname,gname);
      dtk = ak*dk;
      xk_new = xk + dtk;
      fk_new = feval(fname,xk_new);
      dfk = abs(fk - fk_new);
      err = max(dfk,norm(dtk));
      k = k + 1;
end
format long
disp('Solution point:')
xs = xk_new
disp('Value of objective function at the solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations required:')
k