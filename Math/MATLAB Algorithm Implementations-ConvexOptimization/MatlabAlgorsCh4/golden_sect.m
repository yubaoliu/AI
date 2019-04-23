%  Program: golden_sect.m
%  Title: Golden-Section Search
%  Description: Implements the golden-section search.
%  Theory: See Practical Optimization Sec. 4.4
%  Input:
%       fname: objective function
%   [xL1,xU1]: initial interval of uncertainty 
%         rou: final range of uncertainty
%  Output:
%    xs: minimum point
%    fs: minimum value of objective function
%     k: number of iterations at convergence
%  Example: 
%  Find the minimum of 
%     fs=sin(xs)
%  over the range 0.5*pi<= xs <= 2*pi.
%  Solution: 
%  Execute the command
%    [xs,fs,k] = golden_sect('sin',0.5*pi,2*pi,1e-6);
%  Notes:
%  1. The program can be applied to any customized function by
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for typical examples). 
%=============================================================
function [xs,fs,k] = golden_sect(fname,xL1,xU1,rou)
disp(' ')
disp('Program golden.sect.m')
% Generate xa1, xb1, f(xa1), and f(xb1)
K = 0.5*(1+sqrt(5));
k = 1;
I1 = xU1 - xL1;
I2 = I1/K;
xak = xU1 - I2;
xbk = xL1 + I2;
fak = feval(fname,xak);
fbk = feval(fname,xbk);
Ik1 = I2;
xLk = xL1;
xUk = xU1;
% Perform Golden-Section search
while Ik1 >= rou & xak <= xbk,
   Ik2 = Ik1/K;
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
% Display results
format long
disp('Minimum point:')
xs = xw
disp('Minimum value of objective function:')
fs = feval(fname,xw)
format short
disp('Number of iterations performed:')
k