% Program: f_ex15_4.m
% Description: Evaluates the objective and constraint functions 
% for Example 15.4.  It is used to test function sqp_general. 
%=======================================
function z = f_ex15_4(x)
x1 = x(1);
x2 = x(2);
ak = x1^2 + x2^2 - 9;
ck = [x1-1; -x1+5; x2-2; -x2+4];
fk = x1^2 + x2;
z = [fk; ak; ck];