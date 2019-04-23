%  Program: fibonacci.m
%  Title: Fibonacci Search
%  Description: Implements the Fibonacci search.
%  Theory: See Practical Optimization Sec. 4.3
%  Input:
%       fname: objective function
%   [xL1,xU1]: initial interval of uncertainty 
%           n: integer that determines the final range of uncertainty
%  Output:
%     xs: minimum point 
%     fs: minimum value of objective function
%      k: number of iterations at convergence
%  Example: 
%  Find the minimum of 
%     f = sin(x)
%  over the range 0.5*pi <= x <= 2*pi.
%  Solution: 
%  Execute the command
%     [xs,fs,k] = fibonacci('sin',0.5*pi,2*pi,30);
%  Notes:
%  1. The command [X,Y,Z,...] = EVAL(s) returns output arguments 
%     from the expression in string s.
%  2. Command 
%       function [xs,fs,k] = fibonacci(fname,xL1,xU1,n)
%     adds a new function fibonacci to MATLAB's vocabulary.
%  3. The program can be applied to any customized function by
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for a typical example).
%========================================================================
function [xs,fs,k] = fibonacci(fname,xL1,xU1,n)
disp(' ')
disp('Program fibonacci.m')
% Generate a Fibonacci sequence
f(1) = 1;
f(2) = 1;
for i = 1:n-1;
   f(i+2) = f(i) + f(i+1);
end
F = f(2:n+1);
% Generate xa1, xb1, f(xa1), and f(xb1)
k = 1;
I1 = xU1-xL1;
I2 = (F(n-1)/F(n))*I1;
xak = xU1 - I2;
xbk = xL1 + I2;
fak = feval(fname,xak);
fbk = feval(fname,xbk);
Ik1 = I2;
xLk = xL1;
xUk = xU1;
% Perform Fibonacci search
while k < n-2 & xak <= xbk,
   Ik2 = (F(n-k-1)/F(n-k))*Ik1;
   if fak >= fbk,
      xLk1 = xak;
      xUk1 = xUk;
      xak1 = xbk;
      xbk1 = xLk1 + Ik2;
      fak1 = fbk;
      fbk1 = feval(fname,xbk1);
      xw = 0.5*(xak1 + xUk1);
   else
      xLk1 = xLk;
      xUk1 = xbk;
      xak1 = xUk1 - Ik2;
      xbk1 = xak;
      fbk1 = fak;
      fak1 = feval(fname,xak1);
      xw = 0.5*(xLk1 + xbk1);
   end
   xLk = xLk1;
   xUk = xUk1;
   xak = xak1;
   xbk = xbk1;
   fak = fak1;
   fbk = fbk1;
   Ik1 = Ik2;
   k = k + 1;
end
xs = 0.5*(xLk + xUk);
fs = feval(fname,xs);
% Display results
format long
disp('Minimum point:')
xs
disp('Minimum value of objective function:')
fs
format short
disp('Number of iterations performed:')
k