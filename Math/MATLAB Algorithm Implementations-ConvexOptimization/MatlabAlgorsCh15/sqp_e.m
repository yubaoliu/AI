% Program: sqp_e.m
% Title: SQP algorithm for nonlinear problems with 
% equality constraints.
% Description: Implements the SQP algorithm (Algorithm 15.1). 
% Theory: See Practical Optimization Sec. 15.2.1.
% Input:
%   faname: function that evaluates the objective 
%           and constraint functions
%    gname: function that evaluates the gradient of the 
%           objective and constraint functions
%    wname: function that evaluates the Hessian of the objective 
%           and constraint functions
%       x0: initial point
%     lmd0: initial Lagrange multiplier
%     epsi: termination tolerance
% Output:
%       xs: solution point
%       fs: value of objective function at xs.
%        k: number of iterations at convergence
% Example:
% Apply Algorithm 15.1 to solve the minimization problem 
% in Example 15.1
% Solution:
% Execute the following commands:
% x0 = [3 1.5 3]'
% lmd0 = [-1 -1]'
% epsi = 1e-8
% [xs,fs,k] = sqp_e('f_ex15_1','g_ex15_1','w_ex15_1',x0,lmd0,epsi)
% =========================================================
function [xs,fs,k] = sqp_e(faname,gname,wname,x0,lmd0,epsi)
disp(' ')
disp('Program sqp_e.m')
xk = x0(:);
lmdk = lmd0;
n = length(xk);
p = length(lmdk);
p1 = p + 1;
k = 0;
d = 1;
ze = zeros(p,p);
while d >= epsi,
    fak = feval(faname,xk);
    Gk = feval(gname,xk);
    ak = fak(2:p1);
    gk = Gk(:,1);
    Ak = Gk(:,2:p1)';
    xlmdk = [xk; lmdk];
    Wk = feval(wname,xlmdk);
    xz = inv([Wk -Ak'; -Ak ze])*[Ak'*lmdk-gk; ak];
    d_x = xz(1:n);
    xk = xk + d_x;
    lmdk = pinv(Ak')*(Wk*d_x+gk);
    k = k + 1;
    d = norm(d_x);
end
format long
disp('Solution point:')
xs = xk
disp('Objective function at the solution point:')
fak = feval(faname,xs);
fs = fak(1)
format short
disp('Number of iterations at convergence:')
k