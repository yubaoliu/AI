% Program: lsearch_han.m
% Description: Implements Han's line search method.
% It is used in function sqp_ie_h. 
% Theory: See Sec. 15.3.1.
%======================================================
function a = lsearch_han(fcname,xk,d_x,dlt,r)
aa = (0:0.01:1)*dlt;
ps = zeros(101,1);
for i = 1:101,
    ai = aa(i);
    xdi = xk + ai*d_x;
    fci = feval(fcname,xdi);
    fi = fci(1);
    cki = fci(2:end);
    ck0 = max(0,-cki);
    ps(i) = fi + r*sum(ck0);
end
[psm,ind] = min(ps);
a = aa(ind);