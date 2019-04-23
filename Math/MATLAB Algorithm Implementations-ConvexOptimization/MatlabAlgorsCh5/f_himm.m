% Program: f_himm.m
% Description: Evaluates the Himmelblau function
% f(x) =(x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
% To Run: Execute the commands:
% x = [x1 x2]' % e.g., x = [1 2]'
% f = f_himm(x)
% ===============================
function f = f_himm(x)
x1 = x(1);
x2 = x(2);
f = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2;