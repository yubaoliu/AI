% Program: least_pth_m.m
% Title: Modified least-pth minimax algorithm.
% Description: Implements the modified least-pth minimax 
% algorithm (Algorithm 8.3) for the design of 
% stable lowpass, highpass, bandpass, and bandstop IIR digital 
% filters.
% Theory: See Practical Optimization Sec. 8.4.
% This MATLAB function requires the following MATLAB functions:
%   bfgs_least_pth.m; eval_f.m; f_least_pth.m; g_least_pth.m;
%   inex_lsearch.m; initialize.m, stabilize.m
% Input:
%   typ: type of filter to be designed. There are four possible 
%      choices:
%      typ = 'lowpass' for lowpass filters
%      typ = 'highpass' for highpass filters
%      typ = 'bandpass' for bandpass filters
%      typ = 'bandstop' for bandstop filters
%   f: passband and stopband edges. There are four choices for 
%      the correspoding four type of filters:
%      f = [fp fa] for passband filters
%      f = [fa fp] for stopband filters
%      f = [fa1 fp1 fp2 fa2] for bandpass filters
%      f = [fp1 fa1 fa2 fp2] for bandstop filters
%   N: order of the IIR filter
%   K: number of frequencies used.
%   l: parameter l in Eq. (8.20)
%   epsi1: termination tolerance for Algorithm 8.3
%   epsi2: termination tolerance for the sub-optimization
%        to be carried out in Step 3 of the algorithm
%   x_ini is optional user-supplied initial point of the form
%           x_ini = [a0 a1 ... aN b1 b2 ... bN]
%        where a0 a1 ... aN are the N+1 numerator coefficients
%        and b1 b2 ... bN are the N denominator coefficients
%        of the transfer function as defined in Eq.(8.2).       
%        If x_ini is not initialized by the user, the program
%        uses function initialize.m to design an equivalent 
%        lowpass, highpass, bandpass, or bandstop FIR filter  
%        and then uses the coefficients of the filter so obtained  
%        as the initial numerator transfer function coefficients
%        a0 a1 ... aN.  The denomenator coefficients, on the 
%        other hand, are initialized as 
%          b = [1 b1 b2 ... bN]
%        where b1=b2=...=bN=0 which means that the poles of the        
%        transfer function are initially placed at the orgin
%        of the z plane.
%        Another option available to the user is to write 
%        additional instructions at the end of initialize.m
%        which would overwrite the present default initialization
%        by some other scheme.
% Output:
%   a: numerator coefficients of the IIR filter obtained
%   b: denominator coefficients of the IIR filter
%   as: numerator coefficients of the stabilized IIR filter
%   bs: denominator coefficients of the stabilized IIR filter
%   ep_max: maximum approximation error in passband(s)
%   ea_max: maximum approximation error in stopband(s)
%   k: number of iterations at convergence
% Example:
% Applying Algorithm 8.3, design a 6th-order lowpass stable
% IIR filter with normalized passband edge = 0.4*pi rad/s
% and stopband edge = 0.6*pi rad/s.
% Solution:
% Execute the commands
%   typ = 'lowpass'
%   f = [0.4 0.6]
%   N = 6
%   K = 70
%   l = 5
%   epsi1 = 1e-6
%   epsi2 = 1e-7
%   [a,b,as,bs,ep_max,ea_max,k]=least_pth_m(typ,f,N,K,l,epsi1,epsi2,x_ini)
% or
%   [a,b,as,bs,ep_max,ea_max,k]=least_pth_m(typ,f,N,K,l,epsi1,epsi2)
% ========================================================================
function[a,b,as,bs,ep_max,ea_max,k]= ...
   least_pth_m(typ,f,N,K,l,epsi1,epsi2,x_ini)
disp(' ')
disp('Program least_pth_m.m')
warning off
L = l;
v1 = strcmp(typ,'lowpass');
v2 = strcmp(typ,'highpass');
v3 = strcmp(typ,'bandpass');
v4 = strcmp(typ,'bandstop');
v = [v1 v2 v3 v4];
if v1 == 1,
   fp = f(1);
   fa = f(2);
   if nargin > 7,
      x_ini = x_ini(:);
      a = x_ini(1:(N+1));
      b = [1; x_ini((N+2):end)];
   else
      fc = 0.5*(fp+fa);
      [x_ini,a,b] = initialize(v,fc,N);
   end
   ffs = fp + 1 - fa;
   Kp = round(K*(fp/ffs));
   Ka = K - Kp;
   M0 = [ones(Kp,1); zeros(Ka,1)];
   Kpmax = round(1024*(fp/ffs));
   Kamax = 1024 - Kp;
   M0max = [ones(Kpmax,1); zeros(Kamax,1)];
   Lp = (2*Kp-6)*L + 2;
   La = (2*Ka-6)*L + 2;
   fp = fp*pi;
   fa = fa*pi;
   opmax = 0:fp/(Kpmax-1):fp;
   oamax = fa:(pi-fa)/(Kamax-1):pi;
   omimax = [opmax(:); oamax(:)];
   dfp = fp/(Lp-1);
   dfa = (pi-fa)/(La-1);
   sp0 = [0 ((0:1:(L-1))+0.5)*dfp];
   sp1 = [sp0 ((L:1:((2*Kp-7)*L-1))+0.5)*dfp];
   sp = [sp1 (((2*Kp-7)*L):1:((2*Kp-6)*L-1)+0.5)*dfp fp];
   sa0 = [fa fa+((0:1:(L-1))+0.5)*dfa];
   sa1 = [sa0 fa+((L:1:((2*Ka-7)*L-1))+0.5)*dfa];
   sa = [sa1 fa+(((2*Ka-7)*L):1:((2*Ka-6)*L-1)+0.5)*dfa pi];
   num = [Kp Ka L Lp La];
   data1 = [sp sa];
   M = [ones(1,Lp) zeros(1,La)];
elseif v2 == 1,
   fa = f(1);
   fp = f(2);
   if nargin > 7,
      x_ini = x_ini(:);
      a = x_ini(1:(N+1));
      b = [1; x_ini((N+2):end)];
   else
      fc = 0.5*(fp+fa);
      [x_ini,a,b] = initialize(v,fc,N);
   end
   ffs = fa + 1 - fp;
   Ka = round(K*(fa/ffs));
   Kp = K - Ka;
   M0 = [zeros(Ka,1); ones(Kp,1)];
   Kamax = round(1024*(fa/ffs));
   Kpmax = 1024 - Kamax;
   M0max = [zeros(Kamax,1); ones(Kpmax,1)];
   La = (2*Ka-6)*L + 2;
   Lp = (2*Kp-6)*L + 2;
   fa = fa*pi;
   fp = fp*pi;
   oamax = 0:fa/(Kamax-1):fa;
   opmax = fp:(pi-fp)/(Kpmax-1):pi;
   omimax = [oamax(:); opmax(:)];
   dfa = fa/(La-1);
   dfp = (pi-fp)/(Lp-1);
   sa0 = [0 ((0:1:(L-1))+0.5)*dfa];
   sa1 = [sa0 ((L:1:((2*Ka-7)*L-1))+0.5)*dfa];
   sa = [sa1 (((2*Ka-7)*L):1:((2*Ka-6)*L-1)+0.5)*dfa fa];
   sp0 = [fp fp+((0:1:(L-1))+0.5)*dfp];
   sp1 = [sp0 fp+((L:1:((2*Kp-7)*L-1))+0.5)*dfp];
   sp = [sp1 fp+(((2*Kp-7)*L):1:((2*Kp-6)*L-1)+0.5)*dfp pi];
   num = [Ka Kp L La Lp];
   data1 = [sa sp];
   M = [zeros(1,La) ones(1,Lp)];
elseif v3 == 1,
   fa1 = f(1);
   fp1 = f(2);
   fp2 = f(3);
   fa2 = f(4);
   if nargin > 7,
      x_ini = x_ini(:);
      a = x_ini(1:(N+1));
      b = [1; x_ini((N+2):end)];
   else
      fc = [0.5*(fa1+fp1) 0.5*(fa2+fp2)];
      [x_ini,a,b] = initialize(v,fc,N);
   end
   ffs = fa1 + fp2 - fp1 + 1 - fa2;
   Ka1 = round(K*(fa1/ffs));
   Kp = round(K*(fp2-fp1)/ffs);
   Ka2 = K - Ka1 - Kp;
   M0 = [zeros(Ka1,1); ones(Kp,1); zeros(Ka2,1)];
   Ka1max = round(1024*(fa1/ffs));
   Kpmax = round(1024*(fp2-fp1)/ffs);
   Ka2max = 1024 - Ka1max - Kpmax;
   M0max = [zeros(Ka1max,1); ones(Kpmax,1); zeros(Ka2max,1)];
   La1 = (2*Ka1-6)*L + 2;
   Lp = (2*Kp-6)*L + 2;
   La2 = (2*Ka2-6)*L + 2;
   fa1 = fa1*pi;
   fp1 = fp1*pi;
   fp2 = fp2*pi;
   fa2 = fa2*pi;
   oa1max = 0:fa1/(Ka1max-1):fa1;
   opmax = fp1:(fp2-fp1)/(Kpmax-1):fp2;
   oa2max = fa2:(pi-fa2)/(Ka2max-1):pi;
   omimax = [oa1max(:); opmax(:); oa2max(:)];
   dfa1 = fa1/(La1-1);
   dfp = (fp2-fp1)/(Lp-1);
   dfa2 = (pi-fa2)/(La2-1);
   wa0 = [0 ((0:1:(L-1))+0.5)*dfa1];
   wa1 = [wa0 ((L:1:((2*Ka1-7)*L-1))+0.5)*dfa1];
   sa1 = [wa1 (((2*Ka1-7)*L):1:((2*Ka1-6)*L-1)+0.5)*dfa1 fa1];
   wp0 = [fp1 fp1+((0:1:(L-1))+0.5)*dfp];
   wp1 = [wp0 fp1+((L:1:((2*Kp-7)*L-1))+0.5)*dfp];
   sp = [wp1 fp1+(((2*Kp-7)*L):1:((2*Kp-6)*L-1)+0.5)*dfp fp2];
   wa0 = [fa2 fa2+((0:1:(L-1))+0.5)*dfa2];
   wa1 = [wa0 fa2+((L:1:((2*Ka2-7)*L-1))+0.5)*dfa2];
   sa2 = [wa1 fa2+(((2*Ka2-7)*L):1:((2*Ka2-6)*L-1)+0.5)*dfa2 pi];
   num = [Ka1 Kp Ka2 L La1 Lp La2];
   data1 = [sa1 sp sa2];
   M = [zeros(1,La1) ones(1,Lp) zeros(1,La2)];
elseif v4 == 1,
   fp1 = f(1);
   fa1 = f(2);
   fa2 = f(3);
   fp2 = f(4);
   if nargin > 7,
      x_ini = x_ini(:);
      a = x_ini(1:(N+1));
      b = [1; x_ini((N+2):end)];
   else
      fc = [fa1 fa2];
      [x_ini,a,b] = initialize(v,fc,N);
   end
   ffs = fp1 + fa2 - fa1 + 1 - fp2;
   Kp1 = round(K*(fp1/ffs));
   Ka = round(K*(fa2-fa1)/ffs);
   Kp2 = K - Kp1 - Ka;
   M0 = [ones(Kp1,1); zeros(Ka,1); ones(Kp2,1)];
   Kp1max = round(1024*(fp1/ffs));
   Kamax = round(1024*(fa2-fa1)/ffs);
   Kp2max = 1024 - Kp1max - Kamax;
   M0max = [ones(Kp1max,1); zeros(Kamax,1); ones(Kp2max,1)];
   Lp1 = (2*Kp1-6)*L + 2;
   La = (2*Ka-6)*L + 2;
   Lp2 = (2*Kp2-6)*L + 2;
   fp1 = fp1*pi;
   fa1 = fa1*pi;
   fa2 = fa2*pi;
   fp2 = fp2*pi;
   op1max = 0:fp1/(Kp1max-1):fp1;
   oamax = fa1:(fa2-fa1)/(Kamax-1):fa2;
   op2max = fp2:(pi-fp2)/(Kp2max-1):pi;
   omimax = [op1max(:); oamax(:); op2max(:)];
   dfp1 = fp1/(Lp1-1);
   dfa = (fa2-fa1)/(La-1);
   dfp2 = (pi-fp2)/(Lp2-1);
   wp0 = [0 ((0:1:(L-1))+0.5)*dfp1];
   wp1 = [wp0 ((L:1:((2*Kp1-7)*L-1))+0.5)*dfp1];
   sp1 = [wp1 (((2*Kp1-7)*L):1:((2*Kp1-6)*L-1)+0.5)*dfp1 fp1];
   wa0 = [fa1 fa1+((0:1:(L-1))+0.5)*dfa];
   wa1 = [wa0 fa1+((L:1:((2*Ka-7)*L-1))+0.5)*dfa];
   sa = [wa1 fa1+(((2*Ka-7)*L):1:((2*Ka-6)*L-1)+0.5)*dfa fa2];
   wp0 = [fp2 fp2+((0:1:(L-1))+0.5)*dfp2];
   wp1 = [wp0 fp2+((L:1:((2*Kp2-7)*L-1))+0.5)*dfp2];
   sp2 = [wp1 fp2+(((2*Kp2-7)*L):1:((2*Kp2-6)*L-1)+0.5)*dfp2 pi];
   num = [Kp1 Ka Kp2 L Lp1 La Lp2];
   data1 = [sp1 sa sp2];
   M = [ones(1,Lp1) zeros(1,La) ones(1,Lp2)];
end
data = [num data1];
H = freqz(a,b,omimax);
emax = max(abs(abs(H) - M0max));
er = 1;
p = 2;
mu = 2;
k = 0;
disp('Initial point:')
x_ini
while er > epsi1 & k < 9,
    omi = eval_f(x_ini,N,v,data);
    x = bfgs_least_pth(x_ini,N,p,K,M0,omi,epsi2);
    a = x(1:(N+1));
    b = [1; x((N+2):end)];
    H = freqz(a,b,omi);
    pw = abs(abs(H) - M0);
    Pk = max(pw);
    H = freqz(a,b,data1);
    es = abs(abs(H) - M);
    emax1 = max(es);
    er1 = abs(emax1 - emax);
    er2 = abs(Pk - emax1);
    er = max(er1,er2);
    kw = k + 1;
    disp(sprintf('Number of iterations completed: %.5g',kw))
    [H,w] = freqz(a,b,1024);
    if v1 == 1,
       np = round(fp*1024/pi);
       na = round((1-fa/pi)*1024);
       ep = abs(abs(H(1:np))-1);
       ep_max = max(ep);
       ea = abs(H((1025-na):1024));
       ea_max = max(ea);
       ew = [ep(:); zeros(1024-np-na,1); ea(:)];
    elseif v2 == 1,
       na = round(fa*1024/pi);
       np = round((1-fp/pi)*1024);
       ea = abs(H(1:na));
       ea_max = max(ea);
       ep = abs(abs(H((1025-np):1024))-1);
       ep_max = max(ep);
       ew = [ea(:); zeros(1024-np-na,1); ep(:)]; 
    elseif v3 == 1,
       na1 = round(1024*fa1/pi);
       nz1 = round(1024*(fp1-fa1)/pi);
       np = round(1024*(fp2-fp1)/pi);
       na2 = round(1024*(pi-fa2)/pi);
       nz2 = 1024 - (na1 + nz1 + np + na2);
       ea1 = abs(H(1:na1));
       ea1_max = max(ea1);
       ep = abs(abs(H((na1+nz1+1):(na1+nz1+np)))-1);
       ep_max = max(ep);
       ea2 = abs(H((1025-na2):1024));
       ea2_max = max(ea2);
       ew = [ea1(:);zeros(nz1,1);ep(:);zeros(nz2,1);ea2(:)]; 
       ea_max = [ea1_max ea2_max];
    elseif v4 == 1,
       np1 = round(1024*fp1/pi);
       nz1 = round(1024*(fa1-fp1)/pi);
       na = round(1024*(fa2-fa1)/pi);
       np2 = round(1024*(pi-fp2)/pi);
       nz2 = 1024 - (np1 + nz1 + na + np2);
       ep1 = abs(abs(H(1:np1))-1);
       ep1_max = max(ep1);
       ea = abs(H((np1+nz1+1):(np1+nz1+na)));
       ea_max = max(ea);
       ep2 = abs(abs(H((1025-np2):1024))-1);
       ep2_max = max(ep2);
       ew = [ep1(:);zeros(nz1,1);ea(:);zeros(nz2,1);ep2(:)];
       ep_max = [ep1_max ep2_max];
    end
    figure(1)
    Ymax=max(ew);
    plot(w,ew)
    axis([0 pi 0 1.2*Ymax])
    xlabel('Frequency, rad/s')
    ylabel('|Error|')
    grid
    title(sprintf('Approximation error after %.5g iterations',kw))
    pause(1)
    figure(2)
    plot(w,20*log10(abs(H)))
    axis([0 pi -100 10])
    xlabel('Frequency, rad/s') 
    ylabel('Gain, dB')
    grid
    title(sprintf('Amplitude response after %.5g iterations',kw))
    pause(1)
    emax = emax1;
    x_ini = x;
    p = p*mu;
    k = k + 1;
end
disp('Largest magnitude of the poles:')
rt_max = max(abs(roots(b)))
if rt_max < 1,
   as = a;
   bs = b;
else
   [as,bs] = stabilize(a,b);
   disp('Largest magnitude of the poles after stabilization:')
   rt_max_new = max(abs(roots(bs)))
end
disp('Maximum approximation error in passband(s):')
ep_max
disp('Maximum approximation error in stopband(s):')
ea_max
format long