%  Program: DSC_search.m
%  Title: Davies-Swann-Campey Search
%  Description: Implements the Davies-Swann-Campey search.
%  Theory: See Practical Optimization Sec. 4.7
%  Input:
%    fname: objective function
%      x01: initial point
%      dt1: initial increment
%        K: scaling constant
%     epsi: optimization tolerance 
%  Output:
%    xs: minimum point
%    fs: minimum of objective function
%     k: number of iterations at convergence
%    ke: number of function evaluations
%  Example: 
%  Find the minimum of 
%     fs=sin(xs)
%  over the range 0.8*pi<= xs <= 2*pi.
%  Solution: 
%  Execute the command
%      [xs,fs,k,ke] = DSC_search('sin',0.8*pi,0.05,0.1,1e-6);
%  Notes:
%  1. The program can be applied to any customized function by
%     defining the function of interest in an m-file (see file 
%     inex_lsearch.m for typical examples). 
%=============================================================
function [xs,fs,k,ke] = DSC_search(fname,x01,dt1,K,epsi)
disp(' ')
disp('Program DSC_search.m')
k = 0;
ke = 0;
dt = dt1;
x0 = x01;
er = dt;
while er >= epsi,
    x1 = x0 + dt;
    x_1 = x0 - dt;
    f0 = feval(fname,x0);
    f1 = feval(fname,x1);
    ke = ke + 2;
    if f0 > f1,
       p = 1;
       n = 1;
       fn_1 = f0;
       fn = f1;
       xn = x1;
       ind = 0;
    else
       f_1 = feval(fname,x_1);
       ke = ke + 1;
       if f_1 < f0,
          p = -1;
          n = 1;
          fn_1 = f0;
          fn = f_1;
          xn = x_1;
          ind = 0;
       else 
          ind = 1;
       end
    end
    if ind == 0,
       while fn <= fn_1,
             n = n + 1;
             fn_2 = fn_1;
             fn_1 = fn;
             xn_1 = xn;
             xn = xn_1 + (2^(n-1))*p*dt;
             fn = feval(fname,xn);
             ke = ke + 1;
       end
       xm = xn_1 - 2^(n-2)*p*dt;
       fm = feval(fname,xm);
       ke = ke + 1;
       if fm >= fn_1,
          x0 = xn_1+(2^(n-2)*p*dt*(fn_2-fm))/(2*(fn_2-2*fn_1+fm));
       else
          x0 = xm+(2^(n-2)*p*dt*(fn_1-fn))/(2*(fn_1-2*fm+fn));
       end
       er = 2^(n-2)*dt;
       dt = K*dt;
    else
       x0 = x0 + dt*(f_1 - f1)/(2*(f_1 - 2*f0 + f1));
       er = dt;
       dt = K*dt;
    end
    k = k + 1;
end
disp('Minimum point:')
format long
xs = x0
disp('Minimum of objective function:')
fs = feval(fname,xs)
format short
disp('Number of iterations performed:')
k
disp('Number of function evaluations:')
ke