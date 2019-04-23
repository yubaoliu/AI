%  Program: gauss_newton.m
%  Title: Gauss-Newton algorithm
%  Description: implements the Gauss-Newton algorithm (Algorithm 5.5). 
%  The line search is carried out by using the inexact line search 
%  algorithm implemented in terms of Algorithm 4.6.
%  Theory: See Practical Optimization Sec. 5.4.
%  Input: 
%   fname: objective function
%   gname: gradient of the objective function
%   jname: Jacobian of the objective function  
%      x0: initial point
%    epsi: optimization tolerance
%  Output:   
%      xs: solution point
%      fs: objective function evaluated at xs
%       k: number of iterations at convergence
%  Example:
%  Find the minimum of the Himmelblau function
%     f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%  with initioal point x0 = [6 6]' and tolerance 
%  epsi = 1e-6.
%  Solution:
%  Execute the command
%    [xs,fs,k] = gauss_newton('f_himm','g_himm','j_himm',[6 6]',1e-6);
%  Notes:
%  1. The program can be applied to any customized function
%     by defining the function of interest and its gradient and
%     Jacobian in m-files.
% ==================================================================
function [xs,fs,k] = gauss_newton(fname,gname,jname,x0,epsi)
disp(' ')
disp('Program gauss_newton.m')
x = x0(:);
n = length(x);
In = eye(n);
k = 1;
xk = x0;
F_k = feval(fname,xk);
gk = feval(gname,xk);
Jk = feval(jname,xk);
Hk = 2*Jk'*Jk + 1e-12*In;
dk = -inv(Hk)*gk;
ak = inex_lsearch(xk,dk,fname,gname);
adk = ak*dk;
er = norm(adk);
while er > epsi,
  xk = xk + adk;
  F_k1 = feval(fname,xk);
  gk = feval(gname,xk);
  Jk = feval(jname,xk);
  Hk = 2*Jk'*Jk + 1e-12*In;
  dk = -inv(Hk)*gk;
  ak = inex_lsearch(xk,dk,fname,gname);
  adk = ak*dk;
  er = abs(F_k1 - F_k);
  k = k + 1;
  F_k = F_k1;
end
format long
disp('Solution point:')
xs = xk + adk
disp('Objective function at the solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k