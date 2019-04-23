% Program: qp_general.m
% Title: Nonfeasible-initialization primal-dual 
% path-following algorithm for convex QP problems 
% with constraints Ae*x = be and Ax >= b.
% Description: Implements a primal-dual path-following
% algorithm that applies the method in Sec. 13.4.3 to 
% the QP problem 
%              minimize 0.5*x'*H*x +x'*p
%              subject to A*x >= b and Ae = be
% Input:     
%         H: positive semidefinite Hessian matrix
%        Ae: full row-rank constraint matrix Ae
%         A: full row-rank constraint matrix A
%  p, be, b: input vectors
%        x0: initial point such that A*x0 - b > 0
%      epsi: tolerance for duality gap
% Output:   
%        xs: solution vector
% Example:
% Find a minimizer of the QP problem
%   minimize 0.5*x'*H*x +x'*p
%   subject to A*x >= b
% where
%   h = [1 -4 2 1]'; H = h*h'; p = [-1 0 7 4]'; 
%   Ae = [1 1 1 1]; be = 4;
%   A = eye(4); b = zeros(4,1);
%   x0 = [1 1 1 1]';
% Solution:
% Execute the commands:
% h = [1 -4 2 1]'; H = h*h'; p = [-1 0 7 4]'; 
% Ae = [1 1 1 1]; be = 4;
% A = eye(4); b = zeros(4,1);
% x0 = [1 1 1 1]';
% epsi = 1e-6;
% xs = qp_general(H,p,Ae,be,A,b,x0,epsi)
% ============================================
function xs = qp_general(H,p,Ae,be,A,b,x0,epsi)
% disp(' ')
% disp('Program qp_general.m')
% Data initialization.
% n = length(x0);
pe = length(be);
q = length(b);
rho_q = q + 1.5*sqrt(q);
a_max = 1-1e-6;
x = x0(:);
y = A*x - b;
lmd = zeros(pe,1);
mu = ones(q,1);
ym = y.*mu;
gap = sum(ym);
k = 0;
% Iteration begins.
while gap > epsi,
  tau = gap/rho_q;
  rd = -H*x - p + Ae'*lmd + A'*mu;
  ra = be - Ae*x;
  ksi = tau - y.*mu;
  ymi = mu./y;
  yksi = ksi./y;
  YMI = diag(ymi);
  G = inv(H + A'*YMI*A);
  Ga = Ae*G*Ae';
  ayk = rd + A'*yksi;
  rt = ra - Ae*G*ayk;
  % Calculate d_x,d_lmd,and d_mu.
  d_lmd = inv(Ga)*rt;
  d_x = G*(ayk + Ae'*d_lmd);
  d_y = A*d_x;
  d_mu = (ksi - (mu.*d_y))./y;
  % Calculate step size.
  ind = find(d_y < 0);
  a_p = min(y(ind)./(-d_y(ind)));
  ind = find(d_mu < 0);
  a_d = min(mu(ind)./(-d_mu(ind)));
  a_k = a_max*min([1 a_p a_d]);
  % Form new iterate. 
  x = x + a_k*d_x;
  mu = mu + a_k*d_mu;
  lmd = lmd + a_k*d_lmd;
  y = A*x - b;
  % Compute duality gap and check convergence.
  ym = y.*mu;
  gap = sum(ym);
  k = k + 1;
end
xs = x;
% fs = 0.5*xs'*(H*xs + 2*p);