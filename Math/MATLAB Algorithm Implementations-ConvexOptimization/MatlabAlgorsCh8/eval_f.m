% Program: eval_f.m
% Description: Evaluates and selects frequencies
% by using the method described in Step 1 of 
% Algorithms 8.3 and 8.4. This MATLAB function is 
% required by MATLAB functions least_pth_m.m and 
% charalambous_m.m.
% Theory: See Practical Optimization Sec. 8.4.
% ===============================================
function omi = eval_f(x,N,v,data)
a = x(1:(N+1));
b = [1; x((N+1):end)];
v1 = v(1);
v2 = v(2);
v3 = v(3);
v4 = v(4);
if v1 == 1,
   num = data(1:5);
   Kp = num(1);
   Ka = num(2);
   L = num(3);
   Lp = num(4);
   La = num(5);
   sp = data(6:(Lp+5));
   sa = data((Lp+6):end);
   H = freqz(a,b,sp);
   omip(1) = sp(1);
   fw = sp(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))-1));
   omip(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Kp-4),
       fw = sp((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))-1));
       omip(i+2) = fw(ind);
   end
   fw = sp((Lp-L):(Lp-1));
   [er,ind] =max(abs(abs(H((Lp-L):(Lp-1)))-1));
   omip(Kp-1) = fw(ind);
   omip(Kp) = sp(Lp);
   H = freqz(a,b,sa);
   omia(1) = sa(1);
   fw = sa(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omia(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Ka-4),
       fw = sa((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))));
       omia(i+2) = fw(ind);
   end
   fw = sa((La-L):(La-1));
   [er,ind] =max(abs(abs(H((La-L):(La-1)))));
   omia(Ka-1) = fw(ind);
   omia(Ka) = sa(La);
   omi = [omip(:); omia(:)];
elseif v2 == 1,
   num = data(1:5);
   Ka = num(1);
   Kp = num(2);
   L = num(3);
   La = num(4);
   Lp = num(5);
   sa = data(6:(La+5));
   sp = data((La+6):end);
   H = freqz(a,b,sa);
   omia(1) = sa(1);
   fw = sa(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omia(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Ka-4),
       fw = sa((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))));
       omia(i+2) = fw(ind);
   end
   fw = sa((La-L):(La-1));
   [er,ind] = max(abs(abs(H((La-L):(La-1)))));
   omia(Ka-1) = fw(ind);
   omia(Ka) = sa(La);
   H = freqz(a,b,sp);
   omip(1) = sp(1);
   fw = sp(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))-1));
   omip(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Kp-4),
       fw = sp((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))-1));
       omip(i+2) = fw(ind);
   end
   fw = sp((Lp-L):(Lp-1));
   [er,ind] =max(abs(abs(H((Lp-L):(Lp-1)))-1));
   omip(Kp-1) = fw(ind);
   omip(Kp) = sp(Lp);
   omi = [omia(:); omip(:)];
elseif v3 == 1,
   num = data(1:7);
   Ka1 = num(1);
   Kp = num(2);
   Ka2 = num(3);
   L = num(4);
   La1 = num(5);
   Lp = num(6);
   La2 = num(7);
   sa1 = data(8:(La1+7));
   sp = data((La1+8):(La1+7+Lp));
   sa2 = data((La1+8+Lp):end);
   H = freqz(a,b,sa1);
   omia1(1) = sa1(1);
   fw = sa1(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omia1(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Ka1-4),
       fw = sa1((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))));
       omia1(i+2) = fw(ind);
   end
   fw = sa1((La1-L):(La1-1));
   [er,ind] =max(abs(abs(H((La1-L):(La1-1)))));
   omia1(Ka1-1) = fw(ind);
   omia1(Ka1) = sa1(end);
   H = freqz(a,b,sp);
   omip(1) = sp(1);
   fw = sp(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))-1));
   omip(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Kp-4),
       fw = sp((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))-1));
       omip(i+2) = fw(ind);
   end
   fw = sp((Lp-L):(Lp-1));
   [er,ind] =max(abs(abs(H((Lp-L):(Lp-1)))-1));
   omip(Kp-1) = fw(ind);
   omip(Kp) = sp(end);
   H = freqz(a,b,sa2);
   omia2(1) = sa2(1);
   fw = sa2(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omia2(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Ka2-4),
       fw = sa2((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))));
       omia2(i+2) = fw(ind);
   end
   fw = sa2((La2-L):(La2-1));
   [er,ind] =max(abs(abs(H((La2-L):(La2-1)))));
   omia2(Ka2-1) = fw(ind);
   omia2(Ka2) = sa2(end);
   omi = [omia1(:); omip(:); omia2(:)];
else
   num = data(1:7);
   Kp1 = num(1);
   Ka = num(2);
   Kp2 = num(3);
   L = num(4);
   Lp1 = num(5);
   La = num(6);
   Lp2 = num(7);
   sp1 = data(8:(Lp1+7));
   sa = data((Lp1+8):(Lp1+7+La));
   sp2 = data((Lp1+8+La):end);
   H = freqz(a,b,sp1);
   omip1(1) = sp1(1);
   fw = sp1(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omip1(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Kp1-4),
       fw = sp1((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))-1));
       omip1(i+2) = fw(ind);
   end
   fw = sp1((Lp1-L):(Lp1-1));
   [er,ind] =max(abs(abs(H((Lp1-L):(Lp1-1)))-1));
   omip1(Kp1-1) = fw(ind);
   omip1(Kp1) = sp1(end);
   H = freqz(a,b,sa);
   omia(1) = sa(1);
   fw = sa(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))));
   omia(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Ka-4),
       fw = sa((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))));
       omia(i+2) = fw(ind);
   end
   fw = sa((La-L):(La-1));
   [er,ind] =max(abs(abs(H((La-L):(La-1)))));
   omia(Ka-1) = fw(ind);
   omia(Ka) = sa(end);
   H = freqz(a,b,sp2);
   omip2(1) = sp2(1);
   fw = sp2(2:(L+1));
   [er,ind] = max(abs(abs(H(2:(L+1)))-1));
   omip2(2) = fw(ind);
   L2 = 2*L;
   for i = 1:(Kp2-4),
       fw = sp2((L+2+(i-1)*L2):(L+1+i*L2));
       [er,ind] = max(abs(abs(H((L+2+(i-1)*L2):(L+1*i*L2)))-1));
       omip2(i+2) = fw(ind);
   end
   fw = sp2((Lp2-L):(Lp2-1));
   [er,ind] =max(abs(abs(H((Lp2-L):(Lp2-1)))-1));
   omip2(Kp2-1) = fw(ind);
   omip2(Kp2) = sp2(end);
   omi = [omip1(:); omia(:); omip2(:)];
end