% Program: g_ex15_6.m
% Description: Evaluates the gradient of the objective and  
% constraint functions for Example 15.6.  It is used to test 
% function sqp_ie_nc.
%=======================================
function z = g_ex15_6(x)
x1 = x(1);
x2 = x(2);
Ak = [-2*x1 1; -1 2*x2];
gk = [2*(x1-2); 2*(x2-1)];
z = [gk Ak'];