% Program: f_ex15_6.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.6.  It is used to test function sqp_ie_nc. 
%=======================================
function z = f_ex15_6(x)
x1 = x(1);
x2 = x(2);
ck = [-x1^2+x2; -x1+x2^2];
fk = (x1-2)^2+(x2-1)^2;
z = [fk; ck];