% Program: powell.m
% Title: Powell's algorithm 
% Description: Implements Powell's algorithm
% described in Algorithm 6.4. 
% Theory: See Practical Optimization Sec. 6.7.
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
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination 
% tolerance epsi = 1e-6.
% Solution:
% Execute the command
%   [xs,fs,k] = powell('f_himm','g_himm',[6 6]',1e-6)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ========================================================
function [xs,fs,k] = powell(fname,gname,x0,epsi)
disp(' ')
disp('Program powell.m')
n = length(x0);
xk = x0;
D = diag(xk);
xkw = xk;
for i = 1:n,
    dkw = D(:,i);
    a0 = inex_lsearch(xkw,dkw,fname,gname);
    a1 = inex_lsearch(xkw,-dkw,fname,gname);
    if a0 >= a1,
        ak = a0;
    else 
        ak = -a1;
    end
    xkw = xkw + ak*dkw;
end
dk = xkw - xk;
a0 = inex_lsearch(xkw,dk,fname,gname);
a1 = inex_lsearch(xkw,-dk,fname,gname);
if a0 >= a1,
   ak = a0;
else 
   ak = -a1;
end
adk = ak*dk;
err = norm(adk);
k = 1;
while  err >= epsi,
    xk = xkw + adk;
    D = [D(:,2:n) dk];
    xkw = xk;
    for i = 1:n,
        dkw = D(:,i);
        a0 = inex_lsearch(xkw,dkw,fname,gname);
        a1 = inex_lsearch(xkw,-dkw,fname,gname);
        if a0 >= a1,
           ak = a0;
        else 
           ak = -a1;
        end
        xkw = xkw + ak*dkw;
     end
     dk = xkw - xk;
     a0 = inex_lsearch(xkw,dk,fname,gname);
     a1 = inex_lsearch(xkw,-dk,fname,gname);
     if a0 >= a1,
        ak = a0;
     else 
        ak = -a1;
     end
     adk = ak*dk;
     err = norm(adk);
     k = k + 1;
end
format long
disp('Solution point:')
xs = xkw + adk
disp('Objective function at solution point:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k