% Program: zangwill_q.m
% Title: Zangwill's algorithm for convex 
% quadratic functions. 
% Description: Implements Zangwill's algorithm
% applied to convex quaratic functions. 
% Theory: See Practical Optimization Sec. 6.7.
% Input:
%   H: matrix H in the objective function
%     f(x) = 0.5*x'*H*x + x'*b + constant
%   b: vector b in the objective function
%   x0: initial point
%   epsi: termination tolerance
%   epsi1: lower bound for the substitution 
%          of d_km.
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%   k: number of iterations at convergence
% Example:
% Find the minimum of the objective function
%   f = 0.5*x'*H*x + x'*b 
% with H = [1 2; 2 5]; b = [1 -1]';
% using initial point x0 = [9 -11]'.
% Solution:
% Execute the commands
%   H = [1 2; 2 5]
%   b = [1 -1]'
%   x0 = [9 -11]'
%   [xs,fs,k] = zangwill_q(H,b,x0,1e-6,0.2)
% Notes:
% 1. The program can be applied to any customized 
% convex quadratic function.
% ==============================================
function [xs,fs,k] = zangwill_q(H,b,x0,epsi,epsi1)
n = length(x0);
disp(' ')
disp('Program zangwill_q.m')
xk = x0;
D = eye(n);
dtk = 1;
xkw = xk;
for i = 1:n,
    dkw = D(:,i);
    al(i) = -((H*xkw+b)'*dkw)/(dkw'*H*dkw);
    xkw = xkw + al(i)*dkw;
end
[akm,m] = max(al);
dk = xkw - xk;
lamk = norm(dk);
ak = -((H*xkw+b)'*dk)/(dk'*H*dk);
adk = ak*dk;
xk = xkw + adk;
err = max(norm(adk),norm(H*xk+b));
k = 1;
while err >= epsi,
    rk = akm*dtk/lamk;
    if  rk > epsi1,
        D(:,m) = dk;
        dtk = rk;
    end
    xkw = xk;
    for i = 1:n,
        dkw = D(:,i);
        al(i) = -((H*xkw+b)'*dkw)/(dkw'*H*dkw);
        xkw = xkw + al(i)*dkw;
     end
     [akm,m] = max(al);
     dk = xkw - xk;
     lamk = norm(dk);
     ak = -((H*xkw+b)'*dk)/(dk'*H*dk);
     adk = ak*dk;
     xk = xkw + adk;
     err = max(norm(adk),norm(H*xk+b));
     k = k + 1;
end
format long
disp('Solution point:')
xs = xk
disp('Objective function at solution point:')
fs = 0.5*(xs'*H*xs) + xs'*b
format short
disp('Number of iterations performed:')
k