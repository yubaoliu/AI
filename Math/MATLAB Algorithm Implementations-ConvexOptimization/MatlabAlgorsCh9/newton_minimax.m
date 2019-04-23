% Program: newton_minimax.m
% Description: Implements the Newton algorithm
% (Algorithm 5.3) for minimizing the objective 
% function in Eq.(9.45).
% This MATKAB function is required for the 
% implementation of MATLAB function fir_minimax.m
% ===============================================
function x = newton_minimax(x0,v,f,p,w,epsi)
if v(1) == 1,
   fa1 = 0;      
   fp1 = 0;
   fp2 = pi*f(1);
   fa2 = pi*f(2);
elseif v(2) == 1,
   fa1 = pi*f(1);
   fp1 = pi*f(2);
   fp2 = pi;
   fa2 = pi;
elseif v(3) == 1,
   fa1 = pi*f(1);
   fp1 = pi*f(2);
   fp2 = pi*f(3);
   fa2 = pi*f(4);
elseif v(4) == 1,
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
p = 2*p;
pr = [p w v fa1 fp1 fp2 fa2];
k = 1;
xk = x0;
n = length(x0);
I = eye(n);
gk = g_minimax(xk,pr);
Hk = h_minimax(xk,pr);
dk = -inv(Hk+0.01*I)*gk;
ak = inex_lsearch(xk,dk,'f_minimax','g_minimax',pr);
adk = ak*dk;
er = norm(adk);
while er >= epsi,
    xk = xk + adk;
    disp(sprintf('Number of iterations completed: %.5g',k))
    gk = g_minimax(xk,pr);
    Hk = h_minimax(xk,pr);
    dk = -inv(Hk+0.01*I)*gk;
    ak = inex_lsearch(xk,dk,'f_minimax','g_minimax',pr);
    adk = ak*dk;
    er = norm(adk);
    h1 = 0.5*xk(2:end);
    h = [flipud(h1); xk(1); h1];
    [F,ff] = freqz(h,1,4096);
    amp = abs(F);
    figure(1)
    subplot(211)
    plot(ff,20*log10(amp))
    axis([0 pi -80 10])
    grid
    title(sprintf('Amplitude response after %.5g iterations',k))
    pause(0.5)
    subplot(212)
    if v(4) == 0,
       ma = max(amp(i1:i2));
       mi = min(amp(i1:i2));
       plot(ff(i1:i2),amp(i1:i2))
       axis([ff(i1) ff(i2) 0.95*mi 1.05*ma])
    else
       ma = max(max(amp(1:i1)),max(amp(i2:end)));
       mi = min(min(amp(1:i1)),min(amp(i2:end)));
       plot(ff(1:i1),amp(1:i1),'-',ff(i2:end),amp(i2:end),'-')
       axis([0 pi 0.95*mi 1.0*ma])
    end
    grid 
    title(sprintf('Passband ripples after %.5g iterations',k))
    pause(0.5)
    k = k + 1;
end
x = xk + adk;
disp(sprintf('Number of iterations completed: %.5g',k))
h1 = 0.5*x(2:end);
h = [flipud(h1); x(1); h1];
[F,ff] = freqz(h,1,4096);
amp = abs(F);
figure(1)
subplot(211)
plot(ff,20*log10(amp))
axis([0 pi -80 10])
grid
title(sprintf('Amplitude response after %.5g iterations',k))
pause(0.5)
subplot(212)
if v(4) == 0,
   ma = max(amp(i1:i2));
   mi = min(amp(i1:i2));
   plot(ff(i1:i2),amp(i1:i2))
   axis([ff(i1) ff(i2) 0.95*mi 1.05*ma])
else
   ma = max(max(amp(1:i1)),max(amp(i2:end)));
   mi = min(min(amp(1:i1)),min(amp(i2:end)));
   plot(ff(1:i1),amp(1:i1),'-',ff(i2:end),amp(i2:end),'-')
   axis([0 pi 0.95*mi 1.05*ma])
end
grid 
title(sprintf('Passband ripples after %.5g iterations',k))
pause(0.5)