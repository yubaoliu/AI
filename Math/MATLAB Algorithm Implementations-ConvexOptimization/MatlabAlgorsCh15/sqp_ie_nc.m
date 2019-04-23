% Program: sqp_ie_nc.m
% Title: Interior-point algorithm for nonconvex problems with 
% inequality constraints.
% Description: Implements the interior-point algorithm for 
% nonconvex problems with inequality constraints (Algorithm 
% 15.6). 
% Theory: See Practical Optimization Section 15.4.3. 
% Input:
%   fcname: function for the evaluation of the objective and 
%           constraint functions
%    gname: function for the evaluation of the gradient of the 
%           objective and constraint functions
%    hname: function for the evaluation of the Hessian of the 
%           Lagrangian defined by Eq. (15.51)
%       x0: initial point for x
%       y0: initial point for y
%     lmd0: initial Lagrange multiplier
%     tau0: initial value of barrier parameter
%       eo: outer tolerance for convergence
%       ei: inner tolerance for convergence
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%    k: number of iterations at convergence
% Example:
% Apply Algorithm 15.6 to solve the minimization problem 
% in Example 15.6.
% Solution:
% Execute the following commands:
% x0 = [-1 2]'
% y0 = [20 20]'
% lmd0 = [1 2]'
% tau0 = 0.001
% ei = 5e-4
% eo = 1e-5
% [xs,fs,k] = sqp_ie_nc('f_ex15_6','g_ex15_6','h_ex15_6',x0,y0,lmd0,tau0,ei,eo)
% ================================================================
function [xs,fs,k] = sqp_ie_nc(fcname,gname,hname,x0,y0,lmd0,tau0,ei,eo)
disp(' ')
disp('Program sqp_ie_nc.m')
% Data initialization.
xk = x0(:);
yk = y0(:);
lmdk = lmd0(:);
tau = tau0;
n = length(x0);
q = length(y0);
q1 = q + 1;
L = 0;
K = 0;
do = 1;
E2 = eye(n);
% Outer-loop begins.
while do >= eo,
   xk0 = xk;
   yk0 = yk;
   lmdk0 = lmdk;
   k = 0;
   di = 1;
   % Inner-loop begins.
   while di >= ei,
     bt = 0;
     fck = feval(fcname,xk);
     ck = fck(2:q1);
     Gk = feval(gname,xk);
     gk = Gk(:,1);
     Ak = Gk(:,2:q1)';
     yki = 1./yk;
     wwk = yki.*lmdk;
     Yki = diag(yki);
     Lk = diag(lmdk);
     xlk = [n;q;xk;lmdk];
     Hk = feval(hname,xlk);
     Nk = Hk + Ak'*(Yki*Lk)*Ak;
     ew = min(eig(Nk));
     if ew <= 0,
        eb = -1.2*ew + 1e-6;
        Nk = Nk + eb*E2;
     end
     Nki = inv(Nk);
     gak = tau*yki - lmdk;
     rk = yk - ck;
     dx = Nki*(-gk + tau*Ak'*yki + Ak'*(wwk.*rk));
     dL = wwk.*(rk - Ak*dx) + gak;
     dy = Ak*dx - rk;
     % Check if the search direction is descent.
     ksi = gk - tau*Ak'*yki;
     sk = -ksi'*Nki*ksi + tau*yki'*rk + ksi'*(Nki*Ak')*(wwk.*rk);
     if sk >= 0,
        bmin  = sk/(norm(rk)^2);
        bt = 10*bmin;
     end
     % Perform a line search.
     awk = max([-dy./yk; -dL./lmdk]);
     ak = 0.95/awk;
     aks = lsearch_merit(fcname,xk,yk,dx,dy,ak,tau,bt);
     xk = xk + aks*dx;
     yk = yk + aks*dy;
     lmdk = lmdk + aks*dL;
     k = k + 1;
     di = aks*(norm(dx) + norm(dy) + norm(dL));
  end
  % Compute barrier parameter tau.
  tsi = q*min(yk.*lmdk)/(yk'*lmdk);
  tau = 0.1*(min([0.05*(1-tsi)/tsi,2]))^3*(yk'*lmdk)/q;
  K = k + K;
  L = L + 1;
  do = norm(xk-xk0) + norm(yk-yk0) + norm(lmdk-lmdk0);
end
xs = xk;
fck = feval(fcname,xk);
fs = fck(1);
k = K;