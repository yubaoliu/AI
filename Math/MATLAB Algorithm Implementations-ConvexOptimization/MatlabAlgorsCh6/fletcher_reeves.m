% Program: fletcher_reeves.m
% Title: Fletcher-Reeves algorithm 
% Description: Implements the Fletcher-Reeves algorithm
% described in Algorithm 6.3. 
% Theory: See Practical Optimization Sec. 6.6.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   x0: initial point
%   epsi: termination tolerance
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Find the minimum of the Himmelblau function
% f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination tolerance
% epsi = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = fletcher_reeves('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ==========================================================
function [xs,fs,k] = fletcher_reeves(fname,gname,x0,epsi)
disp(' ')
disp('Program fletcher_reeves.m')
k = 1;
xk = x0;
gk = feval(gname,xk);
dk = -gk;
g2 = gk'*gk;
ak = inex_lsearch(xk,dk,fname,gname);
adk = ak*dk;
err = norm(adk);
while  err >= epsi,
    xk = xk + adk;
    gk = feval(gname,xk);
    g2_new = gk'*gk;
    bk = g2_new/g2;
    dk = -gk + bk*dk;
    g2 = g2_new;
    ak = inex_lsearch(xk,dk,fname,gname);
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