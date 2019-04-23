% Program: bfgs_mless.m
% Title: Quasi_Newton Memoryless BFGS algorithm
% Description: Implements the quasi-Newton algorithm
% with memoryless Broyden-Fletcher-Goldfarb-Shanno
% updating formula described in Problem 7.16.
% Theory: See Practical Optimization Secs. 7.6 and 7.10,
% and Problem 7.16.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   x0: initial point
%   epsi1: termination tolerance
% Output:
%   xs: solution point
%   fs: value of objective function at point xs.
%   k: number of iterations required
% Example:
% Find the minimum of the Himmelblau function
% f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination 
% tolerance epsi1 = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = bfgs_mless('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ========================================================
function [xs,fs,k] = bfgs_mless(fname,gname,x0,epsi1)
disp(' ')
disp('Program bfgs_mless.m')
k = 1;
xk = x0;
fk = feval(fname,xk);
gk = feval(gname,xk);
dk = -gk;
ak = inex_lsearch(xk,dk,fname,gname);
dtk = ak*dk;
xk_new = xk + dtk;
fk_new = feval(fname,xk_new);
dfk = abs(fk - fk_new);
err = max(dfk,norm(dtk));
while err >= epsi1,
      gk_new = feval(gname,xk_new);
      gmk = gk_new - gk;
      et4 = dtk'*gmk;
      et1 = (dtk'*gk_new)/et4;
      et2 = (gmk'*gk_new)/et4;
      et3 = (1+(gmk'*gmk)/et4)*et1;
      fk = fk_new;
      gk = gk_new;
      xk = xk_new;
      dk = -gk + et1*gmk + (et2-et3)*dtk;
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