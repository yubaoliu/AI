% f_himm_p.m
function z = f_himm_p(x,p)
x1 = x(1);
x2 = x(2);
a = p(1);
b = p(2);
z = (x1^2 + x2 - a^2)^2 + (x1 + x2^2 - b^2)^2;