% Program: g_ex15_4.m
% Description: Evaluates the gradient of the objective and  
% constraint functions for Example 15.4. It is used to test 
% function sqp_general.
%=======================================
function z = g_ex15_4(x)
x1 = x(1);
x2 = x(2);
Aek = [2*x1 2*x2];
Aik = [1 0; -1 0; 0 1; 0 -1];
gk = [2*x1; 1];
z = [gk Aek' Aik'];