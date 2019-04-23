%  Program: quad_inter.m
%  Title: Quadratic Interpolation Search
%  Description: Implements the quadratic interpolation search.
%  Theory: See Practical Optimization Sec. 4.5
%  Input:
%       fname: objective function
%   [xL1,xU1]: initial interval of uncertainty 
%        epsi: optimization tolerance 
%  Output:
%     xs: minimum point
%     fs: minimum value objective function
%      k: number of iterations at convergence
%  Example: 
%  Find the minimum of 
%     fs=sin(xs)
%  over the range 0.5*pi<= xs <= 2*pi.
%  Solution: 
%  Execute the command
%      [xs,fs,k] = quad_inter('sin',0.5*pi,2*pi,1e-6);
%  Notes:
%  1. The program can be applied to any customized function
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for typical examples). 
%==============================================================
function [xs,fs,k] = quad_inter(fname,x1,x3,epsi)
disp(' ')
disp('Program quad_inter.m')
x0bar = 1e99;
% Perform the first iteration
x2 = 0.5*(x1 + x3);
f1 = feval(fname,x1);
f2 = feval(fname,x2);
f3 = feval(fname,x3);
z1 = (x2 - x3)*f1;
z2 = (x3 - x1)*f2;
z3 = (x1 - x2)*f3;
z4 = (x2 + x3)*z1+(x3 + x1)*z2+(x1 + x2)*z3;
xbar = z4/(2*(z1 + z2 + z3));
fbar = feval(fname,xbar);
d = abs(x0bar - xbar);
k = 1;
% Perform quadratic imterpolation based search
while d >= epsi,
   if x1 < xbar & xbar < x2,
      if fbar <= f2,
        x3 = x2;
        f3 = f2;
        x2 = xbar;
        f2 = fbar;
      else
        x1 = xbar;
        f1 = fbar;
      end
    elseif x2 < xbar & xbar < x3,
      if fbar <= f2,
         x1 = x2;
         f1 = f2;
         x2 = xbar;
         f2 = fbar;
      else 
         x3 = xbar;
         f3 = fbar;
      end
   end
   x0bar=xbar;
   z1 = (x2 - x3)*f1;
   z2 = (x3 - x1)*f2;
   z3 = (x1 - x2)*f3;
   z4 = (x2 + x3)*z1 + (x3 + x1)*z2 + (x1 + x2)*z3;
   xbar = z4/(2*(z1 + z2 + z3));
   fbar = feval(fname,xbar);
   d = abs(x0bar - xbar);
   k = k+1;
end
% Display results
format long
disp('Minimum point:')
xs = xbar
disp('Minimum value of objective function:')
fs = feval(fname,xbar)
format short
disp('Number of iterations performed:')
k