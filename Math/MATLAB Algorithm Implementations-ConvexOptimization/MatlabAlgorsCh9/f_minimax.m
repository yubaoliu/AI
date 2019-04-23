% Program: f_minimax.m
% Description: evaluates the objective function 
% f(h) in Eq. (9.45). This MATLAB function is 
% required by MATLAB function fir_minimax.m.
% ==============================================
function f = f_minimax(x,pr)
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
x = x(:);
n = length(x);
N = 1000;
ind = 0:1:n-1;
ind = ind(:);
if v4 == 0,
   L = fa1 + (fp2-fp1) + (pi-fa2);
   vw1 = 0;
   df1 = 0;
   if fa1 ~= 0,
      N1 = ceil((fa1/L)*N);
      df1 = fa1/(N1-1);
      f1 = 0:df1:fa1;
      vw1 = zeros(N1,1);
      for i = 1:N1-1;
          fw = f1(i) + 0.5*df1;
          vw1(i) = (x'*cos(fw*ind))^2;
      end
   end
   vw2 = 0;
   df2 = 0;
   if fp1 ~= fp2,
      N2 = ceil(((fp2-fp1)/L)*N);
      df2 = (fp2-fp1)/(N2-1);
      f2 = fp1:df2:fp2;
      vw2 = zeros(N2,1);
      for i = 1:N2-1;
          fw = f2(i) + 0.5*df2;
          vw2(i) = (x'*cos(fw*ind)-1)^2;
      end
   end
   vw3 = 0;
   df3 = 0; 
   if fa2 ~= pi,
      N3 = ceil(((pi-fa2)/L)*N);
      df3 = (pi-fa2)/(N3-1);
      f3 = fa2:df3:pi;
      vw3 = zeros(N3,1);
      for i = 1:N3-1;
          fw = f3(i) + 0.5*df3;
          vw3(i) = (x'*cos(fw*ind))^2;
      end
   end
   p2 = round(0.5*p);
   mn = mean([vw1; vw2; vw3]);
   z1 = w*df1*((vw1/mn).^p2);
   z2 = df2*((vw2/mn).^p2);
   z3 = w*df3*((vw3/mn).^p2);
   f = sqrt(mn)*(sum(z1) + sum(z2) + sum(z3))^(1/p);
else
   L = fp1 + (fa2-fa1) + (pi-fp2);
   N1 = ceil((fp1/L)*N);
   df1 = fp1/(N1-1);
   f1 = 0:df1:fp1;
   vw1 = zeros(N1,1);
   for i = 1:N1-1;
       fw = f1(i) + 0.5*df1;
       vw1(i) = (x'*cos(fw*ind)-1)^2;
   end
   N2 = ceil(((fa2-fa1)/L)*N);
   df2 = (fa2-fa1)/(N2-1);
   f2 = fa1:df2:fa2;
   vw2 = zeros(N2,1);
   for i = 1:N2-1;
       fw = f2(i) + 0.5*df2;
       vw2(i) = (x'*cos(fw*ind))^2;
   end
   N3 = ceil(((pi-fp2)/L)*N);
   df3 = (pi-fp2)/(N3-1);
   f3 = fp2:df3:pi;
   vw3 = zeros(N3,1);
   for i = 1:N3-1;
       fw = f3(i) + 0.5*df3;
       vw3(i) = (x'*cos(fw*ind)-1)^2;
   end
   p2 = round(0.5*p);
   mn = mean([vw1; vw2; vw3]);
   z1 = df1*((vw1/mn).^p2);
   z2 = w*df2*((vw2/mn).^p2);
   z3 = df3*((vw3/mn).^p2);
   f = sqrt(mn)*(sum(z1) + sum(z2) + sum(z3))^(1/p);
end