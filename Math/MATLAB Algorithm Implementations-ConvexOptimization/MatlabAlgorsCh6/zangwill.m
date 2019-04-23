% Program: zangwill.m
% Title: Zangwill's algorithm 
% Description: Implements Zangwill's algorithm
% described in Algorithm 6.5. 
% Theory: See Practical Optimization Sec. 6.7.
% Input:
%   fname: objective function
%   gname: gradient of the objective function
%   x0: initial point
%   epsi: termination tolerance
%   epsi1: lower bound for the substitution 
%          of d_km.
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Find the minimum of the Himmelblau function
%   f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% using initial point x0 = [6 6]' and termination 
% tolerance epsi = 1e-4.
% Solution:
% Execute the command
%   [xs,fs,k] = zangwill('f_himm','g_himm',[6 6]',1e-4,0.2)
% Notes:
% 1. The program can be applied to any customized function
% by defining the function of interest, and its gradient.
% ========================================================
function [xs,fs,k] = zangwill(fname,gname,x0,epsi,epsi1)
disp(' ')
disp('Program zangwill.m')
n = length(x0);
xk = x0;
D = eye(n);
dtk = 1;
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
    al(i) = ak;
    xkw = xkw + ak*dkw;
end
[akm,m] = max(al);
dk = xkw - xk;
lamk = norm(dk);
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
    rk = akm*dtk/lamk;
    if  rk > epsi1,
        D(:,m) = dk;
        dtk = rk;
    end
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
        al(i) = ak;
        xkw = xkw + ak*dkw;
     end
     [akm,m] = max(al);
     dk = xkw - xk;
     lamk = norm(dk);
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