% Program: fletcher_switch.m
% Title: Quasi_Newton algorithm with the Fletcher switch 
% method
% Description: Implements the quasi-Newton algorithm
% with the Fletcher switch method applied to the 
% Broyden family updating formula.
% Theory: See Practical Optimization Secs. 7.8 and 7.10.
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
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination 
% tolerance epsi1 = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = fletcher_switch('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ==========================================================
function [xs,fs,k] = fletcher_switch(fname,gname,x0,epsi1)
disp(' ')
disp('Program fletcher_switch.m')
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
while err >= epsi1,
      gk_new = feval(gname,xk_new);
      gmk = gk_new - gk;
      D = dtk'*gmk;
      if D <= 0,
         Sk = I;
      else
         sg = Sk*gmk;
         sw0 = gmk'*sg;
         sw1 = (1+sw0/D)/D;
         sw2 = dtk*dtk';
         sw3 = sg*dtk';
         Sk = Sk + sw1*sw2 - (sw3'+sw3)/D;
      end
      if D - gmk'*Sk*gmk > 0,
         Sk = Sk + sw2/D - (sg*sg')/sw0;
      end
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
disp('Number of iterations at convergence:')
k