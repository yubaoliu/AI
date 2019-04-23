% Program: fir_minimax.m
% Title: Minimax design of linear-phase FIR filters.
% Description: Implements the minimax algorithm for the design
% of near-equiripple linear-phase lowpass, highpass, bandpass, 
% and bandstop FIR digital filters.
% This MATLAB function requires the following MATLAB functions:
%   f_minimax.m, g_minimax.m, h_minimax.m, bfgs_minimax.m,
%   inex_lsearch.m, newton_minimax.m
% Theory: See Practical Optimization Sec. 9.4.2.
% Input:
%      N: order of the IIR filter, N must be an even integer.
%    typ: type of filter to be designed. There are four possible 
%         choices:
%         typ = 'lowpass' for lowpass filters
%         typ = 'highpass' for highpass filters
%         typ = 'bandpass' for bandpass filters
%         typ = 'bandstop' for bandstop filters
%      f: passband and stopband edges. There are four choices for 
%         the correspoding four type of filters:
%         f = [fp fa] for passband filters
%         f = [fa fp] for stopband filters
%         f = [fa1 fp1 fp2 fa2] for bandpass filters
%         f = [fp1 fa1 fa2 fp2] for bandstop filters
%      p: power 2p is used in the least-pth minimization
%      w: a scalar weight for stopband(s)
%   epsi: termination tolerance for the minimization
%        of the least-pth objective fumction in Eq. (9.45)
% choice: choice = 1 for BFGS algorithm, and 
%         choice = 2 for Newton algorithm.
% Output:
%   h: impulse response of the Nth-order FIR filter
% Example:
% Design a 36th-order linear-phase lowpass FIR filter 
% with normalized passband edge = 0.4*pi rad/s and
% stopband edge = 0.6*pi rad/s by using the least-
% squares method described in Sec. 9.4.2.
% Solution:
% Execute the commands
%   N = 36
%   typ = 'lowpass'
%   f = [0.45 0.55]
%   p = 30
%   w = 1
%   epsi = 1e-6
%   choice = 2
%   h = fir_minimax(N,typ,f,p,w,epsi,choice)
% ===============================================
function h = fir_minimax(N,typ,f,p,w,epsi,choice)
disp(' ')
disp('Program fir_minimax.m')
h0 = fir_least_squares(N,typ,f,w);
pause(0.5)
n1 = 1 + round(0.5*N);
x0 = [h0(n1);2*h0(n1+1:end)];
v1 = strcmp(typ,'lowpass');
v2 = strcmp(typ,'highpass');
v3 = strcmp(typ,'bandpass');
v4 = strcmp(typ,'bandstop');
v = [v1 v2 v3 v4];
if v1 == 1,
   fa1 = 0;      
   fp1 = 0;
   fp2 = pi*f(1);
   fa2 = pi*f(2);
elseif v2 == 1,
   fa1 = pi*f(1);
   fp1 = pi*f(2);
   fp2 = pi;
   fa2 = pi;
elseif v3 == 1,
   fa1 = pi*f(1);
   fp1 = pi*f(2);
   fp2 = pi*f(3);
   fa2 = pi*f(4);
elseif v4 == 1,
   fp1 = pi*f(1);
   fa1 = pi*f(2);
   fa2 = pi*f(3);
   fp2 = pi*f(4);
end
if fp1 == 0,
   i1 = 1;
else
   i1 = round(4096*fp1/pi);
end
i2 = round(4096*fp2/pi);
if choice == 1,
   x = bfgs_minimax(x0,v,f,p,w,epsi);
elseif choice == 2,
   x = newton_minimax(x0,v,f,p,w,epsi);
end
h1 = 0.5*x(2:end);
h = [flipud(h1); x(1); h1];
subplot(211)
title('Amplitude Response')
xlabel('Frequency, rad/s')
ylabel('Gain, dB')
subplot(212)
title('Passband Ripple')
xlabel('Frequency, rad/s')
ylabel('Gain')
format long