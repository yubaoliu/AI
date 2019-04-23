% Program: lp_path_sf.m
% Title: Primal-dual path-following algorithm 
%    for standard-form LP problems.
% Description: implements Algorithm 12.4 for 
%    standard-form LP problems.
% Theory: See Practical Optimization Sec. 12.5.1.
% Input:  
%         A,c: input data matrices
% x0,lmd0,mu0: strictly feasible initial point
%         rou: parameter involved in determining the value
%              of tau, rou is required to be no less than
%              square root of n (the dimension of vector c)       
%        epsi: tolerance of duality gap.
% Output:        
%          xs: a minimizer
%          fs: objective function at xs
%           k: number of iteration at convergence
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
% Execute the command
% A = [1 1 1]
% c = [-2 1 -3]'
% x0 = [0.5 0.4 0.1]'
% lmd0 = -4
% mu0 = [2 5 1]'
% rou = 2
% espi = 1e-6
% [xs,fs,k] = lp_path_sf(A,c,x0,lmd0,mu0,rou,epsi)
% =====================================================
function [xs,fs,k] = lp_path_sf(A,c,x0,lmd0,mu0,rou,epsi)
disp('    ')
disp('Program lp_path_sf.m')
% Data initialization
n = length(x0);
rou_n = n + rou;
x = x0(:);
lmd = lmd0(:);
mu = mu0(:);
gap = mu'*x;
k = 0;
% Iteration begins
while gap > epsi,
  % Compute Y, y, etc.
  tau = gap/rou_n;
  mui = 1./mu;
  D = diag(x./mu);
  X = diag(x);
  Y = inv(A*D*A');
  y = x - tau*mui;
  yay = Y*(A*y);
  % Calculate d_x,d_lmd, and d_mu
   d_lmd = yay;
   d_mu = -A'*d_lmd;
   d_x = -D*d_mu-y;
  % Calculate step size
  a_p = [];
  a_d = [];
   for i = 1:n,
    if d_x(i) < 0,
     a_p = [a_p x(i)/(-d_x(i))];
    end
    if d_mu(i) < 0,
     a_d=[a_d mu(i)/(-d_mu(i))];
    end
   end
   a_k=(1-1e-6)*min([min(a_p) min(a_d)]);
   % Form new iterate 
    x = x + a_k*d_x;
    lmd = lmd + a_k*d_lmd;
    mu = mu + a_k*d_mu;
   % Compute duality gap and check convergence
   gap = mu'*x;
 % Update iteration index, etc.
 k = k+1;
end
xs = x;
fs = c'*xs;