% Program: sqp_ie_p.m
% Title: Modified SQP algorithm for nonlinear 
% problems with inequality constraints, where
% the line search is carried out by using Powell's
% method.
% Description: Implements the SQP algorithm (Algorithm 
% 15.3) with Powell's line search. 
% Theory: See Practical Optimization Sec. 15.3.2.
% Input:
%   fcname: objective and constraint functions
%    gname: gradient of the objective and constraint
%           functions
%       x0: initial point
%      mu0: initial Lagrange multiplier
%     epsi: termination tolerance
% Output:
%   xs: solution point
%   fs: objective function evaluated at xs.
%    k: number of iterations at convergence
% Example:
% Apply Algorithm 15.3 with Powell's line search to solve 
% the minimization problem in Example 15.2.
% Solution:
% Execute the following commands:
% x0 = [1 0.5 2 3]'
% mu0 = [1 1]'
% epsi = 1e-5
% [xs,fs,k] = sqp_ie_p('f_ex15_2','g_ex15_2',x0,mu0,epsi)
% =======================================================
function [xs,fs,k] = sqp_ie_p(fcname,gname,x0,mu0,epsi)
disp(' ')
disp('Program sqp_ie_p.m')
xk = x0(:);
muk = mu0(:);
n = length(x0);
q = length(muk);
q1 = q + 1;
In = eye(n);
Yk = In;
fck = feval(fcname,xk);
ck = fck(2:q1);
Gk = feval(gname,xk);
gk = Gk(:,1);
Ak = Gk(:,2:q1)';
k = 0;
d = 1;
while d >= epsi,
    % d_x = quadprog(Yk,gk,-Ak,ck);
    d_x = qp_path_ie(Yk,gk,Ak,-ck,zeros(n,1),epsi);
    muk = pinv(Ak')*(Yk*d_x + gk);
    ala = lsearch_powell(fcname,xk,d_x,muk);
    xk = xk + ala*d_x;
    Gk = feval(gname,xk);
    gk1 = Gk(:,1);
    Ak1 = Gk(:,2:q1)';
    gama_k = (gk1-gk) - (Ak1-Ak)'*muk;
    qk = Yk*d_x;
    dg = d_x'*gama_k;
    ww = d_x'*qk;
    if dg >= 0.2*ww,
       thet = 1;
    else
       thet = 0.8*ww/(ww-dg);
    end
    eta = thet*gama_k + (1-thet)*qk;
    phi = 1/ww;
    cta = 1/(d_x'*eta);
    Yk = Yk + cta*(eta*eta') - phi*(qk*qk');
    Ak = Ak1;
    gk = gk1;
    fck = feval(fcname,xk);
    ck = fck(2:q1);
    k = k + 1;
    d = norm(d_x);
end
format long
disp('Solution point:')
xs = xk
disp('Objective function at the solution point:')
fs = fck(1)
format short
disp('Number of iterations at convergence:')
k