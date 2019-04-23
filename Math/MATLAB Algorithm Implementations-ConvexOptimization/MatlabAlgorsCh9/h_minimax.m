% Program: h_minimax.m
% Description: evaluates Hessian of the 
% objective function f(h) using Eqs.(9.49c)
% -(9.49f). 
% This MATLAB function is required by 
% MATLAB function fir_minimax.m.
% =========================================
function H = h_minimax(x,pr)
p = pr(1);
w = pr(2);
v1 = pr(3);
v2 = pr(4);
v3 = pr(5);
v4 = pr(6);
fa1 = pr(7);
fp1 = pr(8);
fp2 = pr(9);
fa2 = pr(10);
p2 = round(0.5*p);
x = x(:);
n = length(x);
N = 1000;
ind = 0:1:n-1;
ind = ind(:);
q1 = zeros(n,1);
Q1 = zeros(n,n);
zw = f_minimax(x,pr);
zwa = zw^2;
zwb = zw^3;
zwc = 1/sqrt(zw);
if v4 == 0,
   L = fa1 + (fp2-fp1) + (pi-fa2);
   if fa1 ~= 0,
      N1 = ceil((fa1/L)*N);
      df1 = fa1/(N1-1);
      f1 = 0:df1:fa1;
      for i = 1:N1-1;
          fw = f1(i) + 0.5*df1;
          c = cos(fw*ind);
          z11 = x'*c;
          z10 = z11^2/zwa;
          z14 = (z10^(p2-2))/zwb;
          z13 = z14*z10*zwa;
          z15 = z11*c;
          q1 = q1 + (zwc*z13*df1)*z15;
          r1 = (z14*df1)*(z15*z15');
          r2 = (z13*df1/(p-2))*(c*c');
          Q1 = Q1 + r1 + r2;
      end
   end
   q2 = zeros(n,1);
   Q2 = zeros(n,n);
   if fp1 ~= fp2,
      N2 = ceil(((fp2-fp1)/L)*N);
      df2 = (fp2-fp1)/(N2-1);
      f2 = fp1:df2:fp2;
      for i = 1:N2-1;
          fw = f2(i) + 0.5*df2;
          c = cos(fw*ind);
          z21 = x'*c - 1;
          z20 = z21^2/zwa;
          z24 = (z20^(p2-2))/zwb;
          z23 = z24*z20*zwa;
          z25 = z21*c;
          q2 = q2 + (zwc*z23*df2)*z25;
          r1 = (z24*df2)*(z25*z25');
          r2 = (z23*df2/(p-2))*(c*c');
          Q2 = Q2 + r1 + r2;
      end
   end
   q3 = zeros(n,1);
   Q3 = zeros(n,n);
   if fa2 ~= pi,
      N3 = ceil(((pi-fa2)/L)*N);
      df3 = (pi-fa2)/(N3-1);
      f3 = fa2:df3:pi;
      for i = 1:N3-1;
          fw = f3(i) + 0.5*df3;
          c = cos(fw*ind);
          z31 = x'*c;
          z30 = z31^2/zwa;
          z34 = (z30^(p2-2))/zwb;
          z33 = z34*z30*zwa;
          z35 = z31*c;
          q3 = q3 + (zwc*z33*df3)*z35;
          r1 = (z34*df3)*(z35*z35');
          r2 = (z33*df3/(p-2))*(c*c');
          Q3 = Q3 + r1 + r2;
      end
   end
   q = w*q1 + q2 + w*q3;
   Qa = (p-2)*(w*Q1 + Q2 + w*Q3);
   Qb = (p-1)*(q*q');
   H = Qa - Qb;
else
   L = fp1 + (fa2-fa1) + (pi-fp2);
   N1 = ceil((fp1/L)*N);
   df1 = fp1/(N1-1);
   f1 = 0:df1:fp1;
   for i = 1:N1-1;
       fw = f1(i) + 0.5*df1;
       c = cos(fw*ind);
       z11 = x'*c -1;
       z10 = z11^2/zwa;
       z14 = (z10^(p2-2))/zwb;
       z13 = z14*z10*zwa;
       z15 = z11*c;
       q1 = q1 + (zwc*z13*df1)*z15;
       r1 = (z14*df1)*(z15*z15');
       r2 = (z13*df1/(p-2))*(c*c');
       Q1 = Q1 + r1 + r2;
   end
   q2 = zeros(n,1);
   Q2 = zeros(n,n);
   N2 = ceil(((fa2-fa1)/L)*N);
   df2 = (fa2-fa1)/(N2-1);
   f2 = fa1:df2:fa2;
   for i = 1:N2-1;
       fw = f2(i) + 0.5*df2;
       c = cos(fw*ind);
       z21 = x'*c;
       z20 = z21^2/zwa;
       z24 = (z20^(p2-2))/zwb;
       z23 = z24*z20*zwa;
       z25 = z21*c;
       q2 = q2 + (zwc*z23*df2)*z25;
       r1 = (z24*df2)*(z25*z25');
       r2 = (z23*df2/(p-2))*(c*c');
       Q2 = Q2 + r1 + r2;
   end
   q3 = zeros(n,1);
   Q3 = zeros(n,n);
   N3 = ceil(((pi-fp2)/L)*N);
   df3 = (pi-fp2)/(N3-1);
   f3 = fp2:df3:pi;
   for i = 1:N3-1;
       fw = f3(i) + 0.5*df3;
       c = cos(fw*ind);
       z31 = x'*c -1;
       z30 = z31^2/zwa;
       z34 = (z30^(p2-2))/zwb;
       z33 = z34*z30*zwa;
       z35 = z31*c;
       q3 = q3 + (zwc*z33*df3)*z35;
       r1 = (z34*df3)*(z35*z35');
       r2 = (z33*df3/(p-2))*(c*c');
       Q3 = Q3 + r1 + r2;
   end
  q = q1 + w*q2 + q3;
  Qa = (p-2)*(Q1 + w*Q2 + Q3);
  Qb = (p-1)*(q*q');
  H = Qa - Qb;
end