% Program: h_ex15_2.m
% Description: Evaluates the Hessian of the objective and 
% constraint functions for Example 15.2.It is used to test 
% function sqp_ie sqp_ie, sqp_ie_c, sqp_ie_h, and sqp_ie_p.    
%=======================================
function z = h_ex15_2(xlk)
lmdk = xlk(7:8);
lmd1 = lmdk(1);
lmd2 = lmdk(2);
h11 = 1+lmd1/2;
h22 = 1+2*lmd1;
h33 = 1+5*lmd2/4;
h34 = 3*lmd2/4;
z = [h11 0 -1 0; 
0 h22 0 -1; 
-1 0 h33 h34; 
0 -1 h34 h33];