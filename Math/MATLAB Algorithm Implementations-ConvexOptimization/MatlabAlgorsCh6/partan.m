% Program: partan.m
% Title: Partan algorithm
% Description: Implements the partan method
% Theory: See Practical Optimization Sec. 6.8.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   x0: initial point
%   epsi: optimization tolerance
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Find the minimum of the Himmelblau function
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using an initial point x0 = [6 6]' and termination
% tolerance epsi = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = partan('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest and its gradient
% in m-files.
% ========================================================
function [xs,fs,k] = partan(fname,gname,x0,epsi)
disp(' ')
disp('Program partan.m')
k = 1;
xk0 = x0;
gk = feval(gname,xk0);
sk = -gk;
ak = inex_lsearch(xk0,sk,fname,gname);
ask = ak*sk;
xk1 = xk0 + ask;
gk = feval(gname,xk1);
sk = -gk;
ak = inex_lsearch(xk1,sk,fname,gname);
ask = ak*sk;
yk1 = xk1 + ask;
dk = yk1 - xk0;
ak = inex_lsearch(xk0,dk,fname,gname);
adk = ak*dk;
xk1_new = xk0 + adk;
xk0 = xk1;
xk1 = xk1_new;
er = norm(adk);
while er >= epsi,
    gk = feval(gname,xk1);
    sk = -gk;
    ak = inex_lsearch(xk1,sk,fname,gname);
    ask = ak*sk;
    yk1 = xk1 + ask;
    dk = yk1 - xk0;
    ak = inex_lsearch(xk0,dk,fname,gname);
    adk = ak*dk;
    xk1_new = xk0 + adk;
    xk0 = xk1;
    xk1 = xk1_new;
    er = norm(adk);
    k = k + 1;
end
format long
disp('Solution point:')
xs = xk1
disp('Objective function at the solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k