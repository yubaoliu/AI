% Program: steep_desc1.m
% Title: Steepest-descent algorithm without line search
% Description: Implements the steepest-descent method
% described in Algorithm 5.1. Parameter alpha is determined 
% using the closed-form formula in Eq. (5.5) instead of
% a line search. 
% Theory: See Practical Optimization Sec. 5.2.4.
% Input: 
%   fname: objective function
%   gname: gradient of the objective function
%   hname: Hessian of the objective function
%      x0: initial point
%    epsi: termination tolerance
% Output:   
%      xs: solution point
%      fs: objective function evaluated at xs.
%       k: number of iterations at convergence
% Example:
% Find the minimum of the Himmelblau function
%    f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination tolerance 
% epsi = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = steep_desc1('f_himm','g_himm','h_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
%    by defining the function of interest, its gradient, and
%    Hessian in m-files.
% ============================================================
function [xs,fs,k] = steep_desc1(fname,gname,hname,x0,epsi)
disp(' ')
disp('Program steep_desc1.m')
k = 1;
xk = x0;
gk = feval(gname,xk);
Hk = feval(hname,xk);
ak = (gk'*gk)/(gk'*Hk*gk);
adk = -ak*gk;
er = norm(adk);
while er >= epsi,
   xk = xk + adk;
   gk = feval(gname,xk);
   Hk = feval(hname,xk);
   ak = (gk'*gk)/(gk'*Hk*gk);
   adk = -ak*gk;
   er = norm(adk);
   k = k + 1;
end
format long
disp('Solution point:')
xs = xk + adk
disp('Objective function at the solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k