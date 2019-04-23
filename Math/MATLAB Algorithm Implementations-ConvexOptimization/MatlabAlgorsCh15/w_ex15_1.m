% Program: w_ex15_1.m
% Description: Evaluates the Hessian of the objective and constraint 
% functions for Example 15.1. It is used to test function sqp_e. 
%=======================================
function W = w_ex15_1(xlmdk)
x = xlmdk(1:3);
lmd = xlmdk(4:5);
x1 = x(1);
x2 = x(2);
x3 = x(3);
lm1 = lmd(1);
lm2 = lmd(2);
w11 = -12*x1^2-2*x2^2-2*x3^2-lm1*12*x1^2-16*lm2;
w12 = -4*x1*x2;
w13 = -4*x1*x3;
w22 = -24*x2^2-2*x1^2-12*lm1*x2^2-28*lm2;
w33 = -12*x3^2-2*x1^2-12*lm1*x3^2-14*lm2;
W = [w11 w12 w13; w12 w22 0; w13 0 w33];