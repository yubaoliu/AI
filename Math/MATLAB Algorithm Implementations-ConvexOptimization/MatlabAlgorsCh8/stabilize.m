% Program: stabilize.m
% Description: Stabilizes an unstable IIR filter 
% by using the method descibed in Example 8.2 
% in terms of Eqs. (8.29) and (8.30). 
% This MATLAB function is required by the following
% MATLAB functions:
%   least_pth.m, charalambous.m, 
%   least_pth_m.m, charalambous_m.m
% Theory: See Practical Optimization Sec. 8.4.
% =================================================
function [as,bs] = stabilize(a,b)
N = length(b) - 1;
rt =roots(b);
bs = 1;
h0 = 1;
for i = 1:N,
    if abs(rt(i)) < 1,
       bs = conv(bs,[1 -rt(i)]);
    else
       bs = conv(bs,[1 -1/rt(i)]);
       h0 = h0*rt(i);
    end
end
H0 = 1/real(h0);
as = a*H0;
bs = real(bs(:));