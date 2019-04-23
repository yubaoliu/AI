% Program: lp_mehrotra.m
% Title: Mehrotra's nonfeasible-initialization 
%    primal-dual path-following algorithm for 
%    standard-form LP problems.
% Description: Implements Algorithm 12.6 for 
%    standard-form LP problems.
% Theory: See Practical Optimization Sec. 12.5.3.
% Input:  
%       A,b,c: input data matrices
% x0,lmd0,mu0: an initial point with x0 > 0 and mu0 > 0
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
% b = 1
% c = [-2 1 -3]'
% x0 = [0.5 0.5 0.5]'
% lmd0 = 0
% mu0 = [1 1 1]'
% espi = 1e-6
% [xs,fs,k] = lp_mehrotra(A,b,c,x0,lmd0,mu0,epsi)
% =====================================================
function [xs,fs,k] = lp_mehrotra(A,b,c,x0,lmd0,mu0,epsi)
disp('    ')
disp('Program lp_mehrotra.m')
% Initialization
n = length(x0);
x = x0(:);
lam = lmd0(:);
mu = mu0(:);
gap = mu'*x;
k = 0;
% Iteration begins
while gap > epsi,
% Compute D, Y, rd using (12.58d,e) and (12.59d,e)
  tau_h = gap/n;
  mui = 1./mu;
  D = diag(x./mu);
  Y = inv(A*D*A');
  rd = c - A'*lam - mu;
% Calculate predictor direction d_x_a,d_lam_a, and 
% d_mu_a using (12.59a-c)
   d_lam_a = Y*(b+A*D*rd);
   d_mu_a = rd - A'*d_lam_a;
   d_x_a = -x - D*d_mu_a;
% Calculate step sizes a_p_aff and a_d_aff using 
% Eq. (12.60)
  a_p = [];
  a_d = [];
   for i = 1:n,
    if d_x_a(i) < 0,
     a_p = [a_p x(i)/(-d_x_a(i))];
    end
    if d_mu_a(i) < 0,
     a_d = [a_d mu(i)/(-d_mu_a(i))];
    end
   end
   a_p_aff = min([1 min(a_p)]);
   a_d_aff = min([1 min(a_d)]);
 % Compute tau_aff using Eq. (12.61)
   tau_aff = ((mu + a_d_aff*d_mu_a)'*(x + a_p_aff*d_x_a))/n;
 % Determine centering parameter sigma_k and tau 
 % using (12.62) and (12.63)
   sigma_k = (tau_aff/tau_h)^3;
   tau = sigma_k*tau_h;
 % Calculate corrector direction d_x_c,d_lam_c, 
 % and d_mu_c using (12.65)
   y = mui.*d_x_a.*d_mu_a - tau*mui;
   d_lam_c = Y*A*y;
   d_mu_c = -A'*d_lam_c;
   d_x_c = -y - D*d_mu_c;
 % Form searh direction d_x, d_lamda, and d_mu 
 % using (12.66)
   d_x = d_x_a + d_x_c;
   d_mu = d_mu_a + d_mu_c;
   d_lamda = d_lam_a + d_lam_c;
 % Determine step sizes a_k_p and a_k_d using (12.68)
   a_p = [];
   a_d = [];
   for i = 1:n,
    if d_x(i) < 0,
     a_p = [a_p x(i)/(-d_x(i))];
    end
    if d_mu(i) < 0,
     a_d = [a_d mu(i)/(-d_mu(i))];
    end
   end
   a_k_p = min([1 0.99*min(a_p)]);
   a_k_d = min([1 0.99*min(a_d)]);
 % Generate next iterate 
    x = x + a_k_p*d_x;
    lam = lam + a_k_d*d_lamda;
    mu = mu + a_k_d*d_mu;
 % Compute duality gap and check convergence
   gap = mu'*x;
 % Update iteration index, etc.
 k = k + 1;
end
xs = x;
fs = c'*x;