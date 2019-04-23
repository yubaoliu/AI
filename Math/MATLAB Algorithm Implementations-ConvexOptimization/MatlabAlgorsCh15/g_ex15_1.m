% Program: g_ex15_1.m
% Description: Evaluates the gradient of the objective and
% constraint functions for Example 15.1. It is used to test 
% function sqp_e. 
%=======================================
function z = g_ex15_1(x)
 x1 = x(1);
 x2 = x(2);
 x3 = x(3);
 gk = [-4*x1^3-2*x1*x2^2-2*x1*x3^2;
 -8*x2^3-2*x1^2*x2;
 -4*x3^3-2*x1^2*x3];
 Ak = [4*x1^3 4*x2^3 4*x3^3; 16*x1 28*x2 14*x3];
 z = [gk Ak'];