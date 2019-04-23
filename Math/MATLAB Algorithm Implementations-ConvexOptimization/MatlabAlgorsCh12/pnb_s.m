% Program: pnb_s.m
% Title: Primal Newton barrier algorithm for standard-form 
%    LP problems.
% Description: Implements Algorithm 12.2 for standard-form 
%    LP problems.
% Theory: See Practical Optimization Sec. 12.4.
% Input:  
%    A, c: input data matrices
%      x0: a strictly feasible initial point
%    tau0: initial value of barrier parameter
%      ei: inner-loop tolerance
%      eo: outer-loop tolerance
% Output:        
%      xs: a minimizer
%      fs: objective function at xs
%       l: number of outer-loop iterations at convergence
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
% tau0 = 0.1
% ei = 1e-3
% eo = 1e-6
% [xs,fs,l] = pnb_s(A,c,x0,tau0,ei,eo)
% =========================================
function [xs,fs,l]= pnb_s(A,c,x0,tau0,ei,eo)
disp('    ')
disp('Program pnb_s.m')
n = length(x0);
x = x0(:);
do = 1;
tau = tau0;
l = 0;
% Outer loop begins
while do > eo,
 di = 1;
 xL = x;
% Inner loop begins
 while di > ei,
  % Calculate search direction
  X2 = diag(xL).^2;
  lamda = inv(A*X2*A')*A*(X2*c-tau*xL);
  d = xL + X2*(A'*lamda-c)/tau;
  % Calculate upper bound alpha_k_bar
  a_w = [];
   for i = 1:n,
    if d(i) < 0,
     a_w = [a_w xL(i)/(-d(i))];
    end
   end
   if length(a_w) == 0,
      a_bar = 1;
   else
      a_bar = 0.95*min(a_w);
   end
   % Line search using the golden-section method
   par = [n;c(:);xL(:);d(:);tau];
   a = golden_sect_m('f_barrier',0,a_bar,par,1e-4);
   % Form new point
   delta = a*d;
   xL0 = xL;
   xL = xL0 + delta;
   % Check (inner-loop) convergence
   di = norm(delta);
 end
 % Check (outer-loop) convergence
 do = norm(x-xL);
 % Update barrier parameter, iteration index, etc.
 l = l + 1;
 x = xL;
 tau = 0.1*tau;
end
xs = x;
fs = c'*xs;