% g_himm_p.m
function z = g_himm_p(x,p)
x1 = x(1);
x2 = x(2);
a = p(1);
b = p(2);
w1 = (x1^2 + x2 - a^2);
w2 = (x1 + x2^2 - b^2);
z1 = 4*w1*x1 + 2*w2;
z2 = 2*w1 + 4*w2*x2;
z = [z1 z2]';