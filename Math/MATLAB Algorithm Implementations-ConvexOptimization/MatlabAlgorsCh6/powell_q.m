% Program: powell_q.m
% Title: Powell's algorithm for convex 
% quadratic functions. 
% Description: Implements Powell's algorithm
% applied to convex quaratic functions. 
% Theory: See Practical Optimization Sec. 6.7.
% Input:
%   H: matrix H in the objective function
%     f(x) = 0.5*x'*H*x + x'*b + constant
%   b: vector b in the objective function
%   x0: initial point
% Output:
%   xs: solution point
% Example:
% Find the minimum of the objective function
%   f = 0.5*x'*H*x + x'*b 
% with H = [1 2; 2 5] and b = [1 -1]' using 
% initial point x0 = [9 -11]'.
% Solution:
% Execute the commands
%   H = [1 2; 2 5]
%   b = [1 -1]'
%   x0 = [9 -11]'
%   xs = powell_q(H,b,x0)
% Notes:
% 1. The program can be applied to any customized 
% convex quadratic function.
% ==============================================
function xs = powell_q(H,b,x0)
disp(' ')
disp('Program powell_q.m')
n = length(x0);
xk = x0;
D = diag(xk);
xkw = xk;
for i = 1:n,
    dkw = D(:,i);
    gkw = H*xkw + b;
    al = -(gkw'*dkw)/(dkw'*H*dkw);
    xkw = xkw + al*dkw;
end
dk = xkw - xk;
gk = H*xkw + b;
ak = -(gk'*dk)/(dk'*H*dk);
adk = ak*dk;
k = 1;
while  k < n,
    xk = xkw + adk;
    D = [D(:,2:n) dk];
    xkw = xk;
    for i = 1:n,
        dkw = D(:,i);
        gkw = H*xkw + b;
        al = -(gkw'*dkw)/(dkw'*H*dkw);
        xkw = xkw + al*dkw;
    end
    dk = xkw - xk;
    gk = H*xkw + b;
    ak = -(gk'*dk)/(dk'*H*dk);
    adk = ak*dk;
    k = k + 1;
end
format long
disp('Solution point:')
xs = xkw + adk
format short