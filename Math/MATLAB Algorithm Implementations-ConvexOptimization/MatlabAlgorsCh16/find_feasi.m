% Program: find_feasi.m
% Description: This function finds a vector 
% x0 = [x1 x2 ... xn]' such that linear inequality 
% constraint A*x0 > b is satisfied. 
% The code was written based on Algorithm 14.3 by applying
% it to the special case where all matrices involved are
% diagonal.
% Theory: See Practical Optimization, Sec. 14.6.2.
% Input:
% A: an m x n real-valued matrix
% b: an m x 1 real-valued vector
% Output: 
% x0: a vector satisfying A*x0 > b
% Example:
% Find a strictly feasible point for constraint Ax > b, 
% where A = [1 0; -1 0; 0 1; -1 -1; -1 -2];
% b = [0 -2 0 -3.5 -6]'.
% Solution:
% Execute the following commands:
% A = [1 0; -1 0; 0 1; -1 -1; -1 -2]
% b = [0 -2 0 -3.5 -6]'
% x0 = find_feasi(A,b)
% ===========================
function x0 = find_feasi(A,b)
[m,n] = size(A);
nz = n + 1;
mz = m + 1;
FFt = [A -b; zeros(1,n) 1];
% Data initialization
e = ones(mz,1);
e1 = ones(1,nz);
In = 1e-4*eye(nz);
xki = e;
Xw = xki*e1;
W = FFt.*Xw;
Fw = W'*W;
qw = FFt'*xki;
x = inv(Fw+In)*qw;
xkp = FFt*x;
vi = min(xkp);
% Iteration begins
while vi <= 1e-4,
   xw = xki.*xkp - e;
   rou = max(abs(xw));
   gk = 1/(1+rou);
   xki = xki - gk*(xw.*xki);
   Xw = xki*e1;
   W = FFt.*Xw;
   Fw = W'*W;
   qw = FFt'*xki;
   x = inv(Fw+In)*qw;
   xkp = FFt*x;
   vi = min(xkp);
end
% Output result
x0 = x(1:n)/x(n+1);