% Program: h_ex15_6.m
% Description: Evaluates the Hessian of the objective and 
% constraint functions for Example 15.6. It is used to test 
% function sqp_ie_nc.
%=======================================
function z = h_ex15_6(xlk)
lmdk = xlk(5:6);
lmd1 = lmdk(1);
lmd2 = lmdk(2);
z = [2+2*lmd1 0; 0 2-2*lmd2];