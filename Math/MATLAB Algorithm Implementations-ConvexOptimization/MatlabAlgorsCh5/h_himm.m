% Program: h_himm.m
% Description: Evaluates the Hessian of the
% Himmelblau function.
% To Run: Execute the commands:
% x = [x1 x2]' % e.g., x = [1 2]'
% H = h_himm(x)
%===============================
function H = h_himm(x)
x1 = x(1);
x2 = x(2);
h11 = 12*x1^2 + 4*x2 - 42;
h12 = 4*(x1 + x2);
h22 = 4*x1 + 12*x2^2 - 26;
H = [h11 h12; h12 h22];