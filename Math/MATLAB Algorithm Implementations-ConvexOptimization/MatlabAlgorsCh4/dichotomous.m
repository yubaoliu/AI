%  Program: dichotomous.m
%  Title: Dichotomous Search
%  Description: Implements the dichotomous search.
%  Theory: See Practical Optimization Sec. 4.2
%  Input:
%       fname: objective function
%   [xL1,xU1]: initial interval of uncertainty 
%         rou: final range of uncertainty
%  Output:
%    xs: minimum point
%    fs: objective function evaluated at xs
%     k: number of iterations at convergence
%  Example: 
%  Find the minimum of 
%     fs=sin(xs)
%  over the range 0.5*pi<= xs <= 2*pi.
%  Solution: 
%  Execute the command
%    [xs,fs,k] = dichotomous('sin',0.5*pi,2*pi,1e-6);
%  Notes:
%  1. The program can be applied to any customized function by
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for typical examples). 
%=============================================================
function [xs,fs,k] = dichotomous(fname,xL1,xU1,rou)
disp(' ')
disp('Program dichotomous.m')
k = 0;
xL = xL1;
xU = xU1;
I = xU - xL;
% Perform dichotomous search
while I > rou,
   x1 = 0.5*(xL + xU);
   epsi = 0.005*I;
   xa = x1 - epsi;
   xb = x1 + epsi;
   fxa = feval(fname,xa);
   fxb = feval(fname,xb);
   if fxa < fxb,
      xU = xb;
   elseif fxa > fxb,
      xL = xa;
   else
      xL = xa;
      xU = xb;
   end
   I = xU - xL;
   k = k + 1;
end
xs = 0.5*(xL + xU);
fs = feval(fname,xs);
% Display results
disp('Minimum point:')
format long
xs
disp('Minimum of objective function:')
fs
format short
disp('Number of iterations:')
k