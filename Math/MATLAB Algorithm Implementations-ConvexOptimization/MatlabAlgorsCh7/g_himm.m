% Program: g_himm.m
% Description: Evaluates the gradient of the
% Himmelblau function
% To Run: Execute the commands:
% % x = [x1 x2]' % e.g., x = [1 2]'
% g = g_himm(x)
% ================================
function g = g_himm(x)
x1 = x(1);
x2 = x(2);
w1 = (x1^2 + x2 - 11);
w2 = (x1 + x2^2 - 7);
g1 = 4*w1*x1 + 2*w2;
g2 = 2*w1 + 4*w2*x2;
g = [g1 g2]';