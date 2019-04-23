% Program: f_barrier.m
% Title: Objective-function evaluation 
% Description: Evaluates objective function 
%   f_tau(x_k + alpha*d_k) given in Eq. (12.26a).
%   The function evaluation is required 
%   in order to perform a line-search in the primal 
%   Newton barrier algorithm for standard-form LP 
%   problems.
% Theory: See Practical Optimization, Sec. 12.4.1.
% Input: 
%   alpha: scalar variable in line search
%     par: parameters involved in the evaluation,
%          which include
%       n: dimention of x
%       c: cost vector of the LP problem
%      xL: present point x
%       d: search direction
%     tau: value of the barrier parameter
% Output:     
%       f: value of f_tau(x_k + alpha*d_k)
% ==============================================
function f = f_barrier(alpha,par)
n = par(1);
c = par(2:n+1);
c = c(:);
xL = par((n+2):(2*n+1));
xL = xL(:);
d = par((2*n+2):(3*n+1));
d = d(:);
tau = par(3*n+2);
zx = xL + alpha*d;
f = c'*zx - tau*sum(log(zx));