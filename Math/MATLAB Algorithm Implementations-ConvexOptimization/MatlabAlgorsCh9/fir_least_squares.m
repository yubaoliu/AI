% Program: fir_least_squares.m
% Title: Least-squares design of linear-phase FIR filters.
% Description: Implements the weighted least-squares algorithm
% for the design of linear-phase lowpass, highpass, bandpass, 
% and bandstop FIR digital filters.
% Theory: See Practical Optimization Sec. 9.4.1.
% Input:
%     N: order of the IIR filter, N must be an even integer.
%   typ: type of filter to be designed. There are four possible 
%        choices:
%        typ = 'lowpass' for lowpass filters
%        typ = 'highpass' for highpass filters
%        typ = 'bandpass' for bandpass filters
%        typ = 'bandstop' for bandstop filters
%    f: passband and stopband edges. There are four choices for 
%       the correspoding four type of filters:
%    f = [fp fa] for passband filters
%    f = [fa fp] for stopband filters
%    f = [fa1 fp1 fp2 fa2] for bandpass filters
%    f = [fp1 fa1 fa2 fp2] for bandstop filters
%    w: a scalar weight for stopband(s)
% Output:
%    h: impulse response of the Nth-order FIR filter
% Example:
% Design a 36th-order linear-phase lowpass FIR filter 
% with normalized passband edge = 0.4*pi rad/s and
% stopband edge = 0.6*pi rad/s by using the least-
% squares method described in Sec. 9.4.1.
% Solution:
% Execute the commands
%   N = 36
%   typ = 'lowpass'
%   f = [0.45 0.55]
%   w = 1
%   h = fir_least_squares(N,typ,f,1)
% =============================================
function h = fir_least_squares(N,typ,f,w)
disp(' ')
disp('Program fir_least_squares.m')
v1 = strcmp(typ,'lowpass');
v2 = strcmp(typ,'highpass');
v3 = strcmp(typ,'bandpass');
v4 = strcmp(typ,'bandstop');
n = N/2;
if v1 == 1,
   a1 = 0;      
   p1 = 0;
   p2 = pi*f(1);
   a2 = pi*f(2);
elseif v2 == 1,
   a1 = pi*f(1);
   p1 = pi*f(2);
   p2 = pi;
   a2 = pi;
elseif v3 == 1,
   a1 = pi*f(1);
   p1 = pi*f(2);
   p2 = pi*f(3);
   a2 = pi*f(4);
elseif v4 == 1,
   p1 = pi*f(1);
   a1 = pi*f(2);
   a2 = pi*f(3);
   p2 = pi*f(4);
end
Q = zeros(n+1,n+1);
b = zeros(n+1,1);
if v4 == 0,
   Q(1,1) = w*(a1-a2+pi)+p2-p1;
   b(1) = p2 - p1;
   for i = 1:n,
       z1 = sin(2*i*a1) - sin(2*i*a2);
       z2 = sin(2*i*p2) - sin(2*i*p1);
       z3 = sin(i*a1) - sin(i*a2);
       z4 = sin(i*p2 )- sin(i*p1);
       Q(i+1,i+1) = 0.5*Q(1,1) + (w*z1 + z2)/(4*i);
       b(i+1) = z4/i;
   end
   for i = 1:(n+1), 
       for j = (i+1):(n+1),
           z5 = sin((i+j-2)*a1)-sin((i+j-2)*a2);
           z6 = sin((i+j-2)*p2)-sin((i+j-2)*p1);
           z7 = sin((i-j)*a1)-sin((i-j)*a2);
           z8 = sin((i-j)*p2)-sin((i-j)*p1);
           Q(i,j) = 0.5*((w*z5+z6)/(i+j-2)+(w*z7+z8)/(i-j));
           Q(j,i) = Q(i,j);
       end
   end
else
   Q(1,1) = p1 + w*(a2-a1) + pi - p2;
   b(1) = p1 + pi - p2;
   for i = 1:n,
       z1 = sin(2*i*p1);
       z2 = sin(2*i*a2) - sin(2*i*a1);
       z3 = -sin(2*i*p2);
       z4 = sin(i*p1) - sin(i*p2);
       Q(i+1,i+1) = 0.5*Q(1,1) + (z1 + w*z2 + z3)/(4*i);
       b(i+1) = z4/i;
   end
   for i = 1:(n+1), 
     for j = (i+1):(n+1),
         z5 = sin((i-j)*p1);
         z6 = w*(sin((i-j)*a2)-sin((i-j)*a1));
         z7 = -sin((i-j)*p2);
         z8 = sin((i+j-2)*p1);
         z9 = w*(sin((i+j-2)*a2)-sin((i+j-2)*a1));
         z10 = -sin((i+j-2)*p2);
         Q(i,j) = (z5+z6+z7)/(2*(i-j)) + (z8+z9+z10)/(2*(i+j-2));
         Q(j,i) = Q(i,j);
     end
   end
end
xs = inv(Q)*b;
h = zeros(N+1,1);
h(n+1) = xs(1);
h(n+2:N+1) = 0.5*xs(2:n+1);
h(1:n) = flipud(h(n+2:N+1));
[H,ff] = freqz(h,1,1024);
amp = abs(H); 
subplot(211)
plot(ff,20*log10(amp))
axis([0 pi -80 10])
xlabel('Frequency, rad/s')
ylabel('Gain, dB')
title('Amplitude response')
grid
subplot(212)
if p1 == 0,
   i1 = 1;
else
   i1 = round(1024*p1/pi);
end
i2 = round(1024*p2/pi);
if v4 == 0,
         ma = max(amp(i1:i2));
         mi = min(amp(i1:i2));
         plot(ff(i1:i2),amp(i1:i2))
         axis([ff(i1) ff(i2) 0.95*mi 1.05*ma])
else
         plot(ff(1:i1),amp(1:i1),'-',ff(i2:end),amp(i2:end),'-')
         ma = max(max(amp(1:i1)),max(amp(i2:end)));
         mi = min(min(amp(1:i1)),min(amp(i2:end)));
         axis([0 pi 0.95*mi 1.05*ma])
end
grid
title('Passband Ripple')
xlabel('Frequency, rad/s')
ylabel('Gain')
format long