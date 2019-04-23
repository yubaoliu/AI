% Program: initialize.m
% Description: Designs an equivalent lowpass, highpass,
% bandpass, or bandstop FIR filter and then uses the coefficients 
% of the filter so obtained as the initial numerator transfer
% function coefficients a0 a1 ... aN.  The denomenator coefficients, 
% on the other hand, are initialized as 
%          b = [1 b1 b2 ... bN]
% where b1=b2=...=bN=0, which means that the poles of the       
% transfer function are initially placed at the orgin of the z 
% plane.
% This MATLAB function is required by the following
% MATLAB functions:
%      least_pth.m, charalambous.m, 
%      least_pth_m.m, charalambous_m.m.
% which can be used to design IIR filters using minimax
% methods. 
% Note: The user can also include additional instructions
% at the end of the file which can overwrite the above  
% default initialization by some other scheme.
% Input
%     v: a four-component vector to specify filter type:
%        v = [1 0 0 0] for lowpass filter
%        v = [0 1 0 0] for highpass filter
%        v = [0 0 1 0] for bandpass filter
%        v = [0 0 0 1] for bandstop filter
%    fc: normalized cutoff frequency: 
%        fc is a scalar between 0 and 1 for lowpass and 
%        highpass filters
%        fc = [fc1 fc2] for bandpass and bandstop filters.
%     N: Order of the IIR filter.
% Output
%    x_ini: a vector of length 2*N + 1, where the first N+1 
%        components are the coefficients of the numerator 
%        a0, a1, ..., aN and the last N components are
%        the coefficients of the deniminator b1, b2, ..., bN.
%     a: a = [a0 a1 ... aN]
%     b: b = [1 b1 b2 ... bN]
%=============================================================
function [x_ini,a,b] = initialize(v,fc,N)
if v(1) == 1,
   a = fir1(N,fc,'low');
elseif v(2) == 1,
   a = fir1(N,fc,'high');
elseif v(3) == 1,
   a = fir1(N,fc,'bandpass');
elseif v(4) == 1,
   a = fir1(N,fc,'stop');
end
a = a(:);
b = [1; zeros(N,1)];
x_ini = [a; zeros(N,1)];
% The user can write additional instructions here to overwrite 
% the above default initialization by some other scheme.
