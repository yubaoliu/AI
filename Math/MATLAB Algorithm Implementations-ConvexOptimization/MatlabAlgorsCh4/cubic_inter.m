%  Program: cubic_inter.m
%  Title: Cubic Interpolation Search
%  Description: Implements the cubic interpolation search.
%  Theory: See Practical Optimization Sec. 4.6
%  Input:
%    fname: objective function
%    gname: gradient of the objective function
%  [x1,x3]: initial interval of uncertainty 
%     epsi: optimization tolerance 
%  Output:
%       xs: minimum point
%       fs: minimum of objective function
%        k: number of iterations at convergence
%       kg: number of 1st-order derivative evaluated
%  Example: 
%  Find the minimum of 
%     fs=sin(xs)
%  over the range 0.8*pi<= xs <= 2*pi.
%  Solution: 
%  Execute the command
%     [xs,fs,k,kg] = cubic_inter('sin','cos',0.8*pi,2*pi,1e-6);
%  Notes:
%  1. The program can be applied to any customized function by
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for typical examples). 
%==============================================================
function [xs,fs,k,kg] = cubic_inter(fname,gname,x1,x3,epsi)
disp(' ')
disp('Program cubic_inter.m')
x0bar = 1e99;
% perform 1st iteration
x2 = 0.5*(x1+x3);
f1 = feval(fname,x1);
f2 = feval(fname,x2);
f3 = feval(fname,x3);
g1 = feval(gname,x1);
kg = 1;
dx12 = x1 - x2;
dx13 = x1 - x3;
sx12 = x1 + x2;
sx13 = x1 + x3;
bt = (f2 - f1 + g1*dx12)/(dx12^2);
gm = (f3 - f1 + g1*dx13)/(dx13^2);
th = (2*x1^2 - x2*sx12)/dx12;
ps = (2*x1^2 - x3*sx13)/dx13;
a3 = (bt - gm)/(th - ps);
a2 = bt - th*a3;
a1 = g1 - 2*a2*x1 - 3*a3*x1^2;
sw = sqrt(a2^2 - 3*a1*a3);
xe1 = (-a2+sw)/(3*a3);
xe2 = (-a2-sw)/(3*a3);
ct = -a2/(3*a3);
if a3 > 0,
    if xe1 > ct,
       xbar = xe1;
    else
       xbar = xe2;
    end
else
    if xe1 < ct,
       xbar = xe1;
    else
       xbar = xe2;
    end
end
fbar = feval(fname,xbar);
d = abs(x0bar-xbar);
k = 1;
% Perform cubic imterpolation based search
while d >= epsi,
%     [fm,ind] = max([f1 f2 f3]);
%     if ind == 1,
%         x1 = xbar;
%         f1 = fbar;
%         g1 = feval(gname,xbar);
%         kg = kg + 1;
%     else
%         x3 = xbar;
%         f3 = fbar;
%     end
   if x1 < xbar & xbar < x2,
      if fbar <= f2,
        x3 = x2;
        f3 = f2;
        x2 = xbar;
        f2 = fbar;
      else
        x1 = xbar;
        f1 = fbar;
        g1 = feval(gname,x1);
        kg = kg + 1;
      end
    elseif x2 < xbar & xbar < x3,
      if fbar <= f2,
         x1 = x2;
         f1 = f2;
         g1 = feval(gname,x1);
         kg = kg + 1;
         x2 = xbar;
         f2 = fbar;
      else 
         x3 = xbar;
         f3 = fbar;
      end
   end
x0bar = xbar;
dx12 = x1 - x2;
dx13 = x1 - x3;
sx12 = x1 + x2;
sx13 = x1 + x3;
bt = (f2 - f1 + g1*dx12)/(dx12^2);
gm = (f3 - f1 + g1*dx13)/(dx13^2);
th = (2*x1^2 - x2*sx12)/dx12;
ps = (2*x1^2 - x3*sx13)/dx13;
a3 = (bt - gm)/(th - ps);
a2 = bt - th*a3;
a1 = g1 - 2*a2*x1 - 3*a3*x1^2;
sw = sqrt(a2^2 - 3*a1*a3);
xe1 = (-a2+sw)/(3*a3);
xe2 = (-a2-sw)/(3*a3);
ct = -a2/(3*a3);
   if a3 > 0,
      if xe1 > ct,
         xbar = xe1;
      else
         xbar = xe2;
      end
   else
      if xe1 < ct,
         xbar = xe1;
      else
         xbar = xe2;
      end
   end
fbar = feval(fname,xbar);
d = abs(x0bar-xbar);
k = k + 1;
end
% Display results
format long
disp('Minimum point:')
xs = real(xbar)
disp('Minimum value of objective function:')
fs = real(fbar)
format short
disp('Number of iterations performed:')
k
disp('Number of 1st-rder derivatives performed:')
kg