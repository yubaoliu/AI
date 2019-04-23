%  Program: j_himm.m
%  Description: Evaluates the Jacobian of 
%  [f1(x) f2(x)]' for the Himmelblau function, 
%  where
%             f1(x) = x1^2 + x2 - 11
%             f2(x) = x1 + x2^2 - 7
%  To Run: Execute the commands:
%    x = [x1 x2]' % e.g., x=[1 2]'
%    J = j_himm(x)
%  ===============================
function J = j_himm(x)
x1 = x(1);
x2 = x(2);
J = [2*x1 1; 1 2*x2];