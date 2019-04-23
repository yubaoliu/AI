% Program: conj_gradient.m
% Title: Conjugate-gradient algorithm 
% Description: Implements the conjugate-gradient
% algorithm described in Algorithm 6.2. 
% Theory: See Practical Optimization Sec. 6.4.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   hname: Hessian of the objective function
%   x0: initial point
%   epsi: termination tolerance
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Find the minimum of the Himmelblau function
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination tolerance
% epsi = 1e-6.
% Solution:
% Execute the command
%   x0 = [6 6]'
%   epsi = 1e-6
% [xs,fs,k]=conj_gradient('f_himm','g_himm','h_himm',x0,epsi)
% Notes:
% 1. In principle, the program can be applied to any 
% customized function by defining the function of interest, 
% and its gradient and Hessian. However, the algorithm may
% fail to work due primarily to the lack of a reliable 
% line-search step. For non-quadratic objective functions,
% the Fletcher-Reeves algorithm is preferred.
% ==========================================================
function [xs,fs,k]=conj_gradient(fname,gname,hname,x0,epsi)
disp(' ')
disp('Program conj_gradient.m')
k = 1;
xk = x0;
gk = feval(gname,xk);
Hk = feval(hname,xk);
dk = -gk;
g2 = gk'*gk;
ak = g2/(dk'*Hk*dk);
adk = ak*dk;
err = norm(adk);
while  err >= epsi,
    xk = xk + adk;
    gk = feval(gname,xk);
    Hk = feval(hname,xk);
    g2_new = gk'*gk;
    bk = g2_new/g2;
    dk = -gk + bk*dk;
    g2 = g2_new;
    ak = g2/(dk'*Hk*dk);
    adk = ak*dk;
    err = norm(adk);
    k = k + 1;
end
format long
disp('Solution point:')
xs = xk + adk
disp('Objective function at solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k