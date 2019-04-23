%  Program: golden_sect_m.m
%  Title: Modified golden-section Search
%  Description: Implements the golden-section search
%     for functions that contain parameters in addition 
%     to design variables.
%  Theory: See Practical Optimization Sec. 4.4
%  Input:
%       fname: objective function
%   [xL1,xU1]: initial interval of uncertainty 
%         par: parameter vector
%         rou: final range of uncertainty
%  Output:
%           x: line search result
%  Example: 
%  Find the minimum of 
%     f(x) = sin(2*a*x + b)
%  over the range 0.5*pi<= x <= 2*pi where
%  a and b are parameters which assume the values
%  a = 0.45 and b = 0.1.
%  Solution: 
%  (1) Prepare a function file:
%      function f = f_sin(x,par)
%      a = par(1);
%      b = par(2);
%      f = sin(2*a*x+b);
%  (2) Execute the command
%      a = 0.45
%      b = 0.1
%      par = [a b]
%      x = golden_sect_m('f_sin',0.5*pi,2*pi,par,1e-6)
%=============================================================
function x = golden_sect_m(fname,xL1,xU1,par,rou)
% Generate xa1, xb1, f(xa1), and f(xb1)
K = 0.5*(1+sqrt(5));
k = 1;
I1 = xU1 - xL1;
I2 = I1/K;
xak = xU1 - I2;
xbk = xL1 + I2;
fak = feval(fname,xak,par);
fbk = feval(fname,xbk,par);
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
      fbk1 = feval(fname,xbk1,par);
      xw = 0.5*(xak1 + xUk1);
   else
      xLk1 = xLk;
      xUk1 = xbk;
      xak1 = xUk1 - Ik2;
      xbk1 = xak;
      fbk1 = fak;
      fak1 = feval(fname,xak1,par);
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
x = xw;