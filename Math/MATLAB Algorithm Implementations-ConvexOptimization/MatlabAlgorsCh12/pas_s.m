% Program: pas_s.m
% Title: Primal affine scaling algorithm for 
%    standard-form LP problems.
% Description: implements Algorithm 12.1 for 
%    standard-form LP problems.
% Theory: See Practical Optimization Sec. 12.3.
% Input:  
%    A, c: input data matrices
%      x0: a strictly feasible initial point
%     gam: constant gamma in formula (12.25a)
%    epsi: tolerance
% Output:        
%      xs: a minimizer
%      fs: objective function at xs
%       k: number of iteration at convergence
% Example:
% Find a minimizer of the standard-form LP problem
%    minimize c'*x 
%    subject to A*x = b, x>= 0
% where
% A = [1 1 1]
% b = 1
% c = [-2 1 -3]'
% The initial point is given by
% x0 = [0.5 0.4 0.1]'
% Solution:
% Execute the commands
% A = [1 1 1]
% c = [-2 1 -3]'
% x0 = [0.5 0.4 0.1]'
% gam = 0.9999
% epsi = 1e-6
% [xs,fs,k] = pas_s(A,c,x0,gam,epsi)
% ========================================
function [xs,fs,k] = pas_s(A,c,x0,gam,epsi)
disp('    ')
disp('Program pas_s.m')
% Data initialization
n = length(x0);
er = 1;
x = x0;
k = 0;
f(1) = c'*x;
while er > epsi,
  % Compute search direction dk
  X2 = diag(x.^2);
  XA = X2*A';
  dk = -(X2-XA*inv(A*XA)*XA')*c;
  % Determine step size alpha_k
  ind = find(dk < 0);
  r = -x(ind)./dk(ind);
  alpha = gam*min(r);
  % Set new iterate, etc.
  x = x + alpha*dk;
  k = k + 1;
  f(k+1) = c'*x;
  % Check convergence
  er = abs(f(k+1)-f(k))/max(1,abs(f(k)));
end
xs = x;
fs = f(end);