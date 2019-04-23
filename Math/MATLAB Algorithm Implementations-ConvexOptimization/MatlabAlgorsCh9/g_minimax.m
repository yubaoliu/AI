% Program: g_minimax.m
% Description: evaluates gradient of the 
% objective function f(h) using Eqs.(9.49a)
% and (9.49b). 
% This MATLAB function is required by 
% MATLAB function fir_minimax.m.
% =========================================
function g = g_minimax(x,pr)
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
L = fa1 + (fp2-fp1) + (pi-fa2);
ind = 0:1:n-1;
ind = ind(:);
zw = f_minimax(x,pr);
zw2 = zw^2;
if v4 == 0,
   z1 = zeros(n,1);
   if fa1 ~= 0,
      N1 = ceil((fa1/L)*N);
      df1 = fa1/(N1-1);
      f1 = 0:df1:fa1;
      for i = 1:N1-1;
          fw = f1(i) + 0.5*df1;
          c = cos(fw*ind);
          z11 = x'*c;
          z13 = (z11^2/zw2)^(p2-1);
          z14 = z11*c; 
          z1 = z1 + (z13*df1)*z14;
      end
      z1 = w*z1;
   end
   z2 = zeros(n,1);
   if fp1 ~= fp2,
      N2 = ceil(((fp2-fp1)/L)*N);
      df2 = (fp2-fp1)/(N2-1);
      f2 = fp1:df2:fp2;
      for i = 1:N2-1;
          fw = f2(i) + 0.5*df2;
          c = cos(fw*ind);
          z21 = x'*c - 1;
          z23 = (z21^2/zw2)^(p2-1);
          z24 = z21*c;      
          z2 = z2 + (z23*df2)*z24;
      end
   end
   z3 = zeros(n,1);
   if fa2 ~= pi,
      N3 = ceil(((pi-fa2)/L)*N);
      df3 = (pi-fa2)/(N3-1);
      f3 = fa2:df3:pi;
      for i = 1:N3-1;
          fw = f3(i) + 0.5*df3;
          c = cos(fw*ind);
          z31 = x'*c;
          z33 = (z31^2/zw2)^(p2-1);
          z34 = z31*c;         
          z3 = z3 + (z33*df3)*z34;
      end
      z3 = w*z3;
   end
   z = z1 + z2 + z3;
   zw1 = 1/zw;
   g = zw1*z(:);
else
   L = fp1 + (fa2-fa1) + (pi-fp2);
   z1 = zeros(n,1);
   N1 = ceil((fp1/L)*N);
   df1 = fp1/(N1-1);
   f1 = 0:df1:fp1;
   for i = 1:N1-1;
       fw = f1(i) + 0.5*df1;
       c = cos(fw*ind);
       z11 = x'*c - 1;
       z13 = (z11^2/zw2)^(p2-1);
       z14 = z11*c; 
       z1 = z1 + (z13*df1)*z14;
   end
   z2 = zeros(n,1);
   N2 = ceil(((fa2-fa1)/L)*N);
   df2 = (fa2-fa1)/(N2-1);
   f2 = fa1:df2:fa2;
   for i = 1:N2-1;
       fw = f2(i) + 0.5*df2;
       c = cos(fw*ind);
       z21 = x'*c;
       z23 = (z21^2/zw2)^(p2-1);
       z24 = z21*c;      
       z2 = z2 + (z23*df2)*z24;
   end
   z2 = w*z2;
   z3 = zeros(n,1);
   N3 = ceil(((pi-fp2)/L)*N);
   df3 = (pi-fp2)/(N3-1);
   f3 = fp2:df3:pi;
   for i = 1:N3-1;
       fw = f3(i) + 0.5*df3;
       c = cos(fw*ind);
       z31 = x'*c - 1;
       z33 = (z31^2/zw2)^(p2-1);
       z34 = z31*c;         
       z3 = z3 + (z33*df3)*z34;
   end
   z = z1 + z2 + z3;
   zw1 = 1/zw;
   g = zw1*z(:);
end