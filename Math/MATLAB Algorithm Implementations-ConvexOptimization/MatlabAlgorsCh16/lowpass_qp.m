% Program: lowpass_qp.m
% Description: This function implements the least-squares
% design of linear-phase lowpass FIR filters described 
% in Sec. 16.2.1.
% Theory: See Practical Optimization, Sec. 16.2.1.
% Input:     
% N: Order of the FIR filter, N must be an even integer.
% p: normalized passband edge between 0 and 1
% a: normalized stopband edge between 0 and 1
% gam: weight for the stopbands
% dp: maximum pass band error
% da: maximum stopband gain
% Mp: number of sample frequencies in passband
% Ma: number of sample frequencies in stopband
% Output:     
% h: impulse response of the bandpass FIR filter.
% Example: 
% Solve the problem in Example 16.1 by using the design
% method described in Sec. 16.2.1.
% Solution: 
% Execute the following command:
% h = lowpass_qp(68,0.43,0.48,1000,0.03,0.026,120,10)
% ==================================================
function h = lowpass_qp(N,p,a,gam,dp,da,Mp,Ma)
disp(' ')
disp('Program lowpass_qp.m')
p = p*pi;
a = a*pi;
n = N/2;
n1 = n + 1;
cc = 0:1:n;
Ql1 = zeros(n1,n1);
Ql2 = Ql1;
bt = zeros(n1,1);
for i = 1:n1,
   for j = 1:n1,
      if i == j,
         if i == 1,
            Ql1(i,i) = p;
            Ql2(i,i) = gam*(pi-a);
         else
            Ql1(i,i) = p/2 + sin(2*(i-1)*p)/(4*(i-1));
            Ql2(i,i) = gam*((pi-a)/2 - sin(2*(i-1)*a)/(4*(i-1)));
         end
      else
         w1 = sin((i-j)*p)/(i-j) + sin((i+j-2)*p)/(i+j-2);
         Ql1(i,j) = 0.5*w1;
         w2 = sin((i-j)*a)/(i-j) + sin((i+j-2)*a)/(i+j-2);
         Ql2(i,j) = -0.5*gam*w2;
      end
   end
end
Q = Ql1 + Ql2;
bt(1) = p;
for i = 2:n1,
   bt(i) = sin((i-1)*p)/(i-1);
end
sp = 0:p/(Mp-1):p;
sa = a:0.1*(pi-a)/(Ma-1):a+0.1*(pi-a);
Ap1 = zeros(Mp,n1);
bp1 = (1+dp)*ones(Mp,1);
bp2 = (-1+dp)*ones(Mp,1);
bp = [bp1; bp2];
for i = 1:Mp,
   Ap1(i,:) = cos(sp(i)*cc);
end
Ap = [Ap1; - Ap1];
Aa1 = zeros(Ma,n1);
ba = da*ones(2*Ma,1);
for i = 1:Ma,
   Ap1(i,:) = cos(sa(i)*cc);
end
Aa = [Aa1; -Aa1];
A = -[Ap; Aa];
b = -[bp; ba];
x0 = find_feasi(A,b);
xls = qp_path_ie(Q,-bt,A,b,x0,1e-6);
hin = xls;
hin(2:n1) = 0.5*xls(2:n1);
h = [flipud(hin(2:n1)); hin];
f = 0:pi/2047:pi;
[H,w] = freqz(h,1,2048);
fp = 0:p/1023:p;
Hp = freqz(h,1,fp);
erp = max(abs((abs(Hp)-1)));
fa = a:(pi-a)/1023:pi;
Ha = freqz(h,1,fa);
era = -max(20*log10(abs(Ha)));
disp('maximum passband ripple is:')
erp
disp('minimum stopband attenuation (in dB) is:')
era
figure(1)
plot(f,20*log10(abs(H)))
axis([0 pi -110 10])
grid
xlabel('(a)  {\it\omega}, rad/s','fontsize',12)
ylabel('Gain, dB','fontsize',12)
figure(2)
plot(fp,abs(Hp)-1)
axis([0 p -0.05 0.05])
grid
xlabel('(b)  {\it\omega}, rad/s','fontsize',12)
ylabel('Gain, dB','fontsize',12)