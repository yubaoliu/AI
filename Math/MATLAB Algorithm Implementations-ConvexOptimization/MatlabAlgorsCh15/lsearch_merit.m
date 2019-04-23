% Program: lsearch_merit.m
% Description: Implements the line search method for merit 
% function in Eq. (15.57).
% It is used in functions sqp_ie_c and sqp_ie_nc.
% Theory: See Sec. 15.4.2.
%======================================================
function a = lsearch_merit(fcname,xk,yk,dx,dy,ak,tau,bt)
aa = (0:0.005:1)*ak;
ps = zeros(101,1);
for i = 1:101,
    ai = aa(i);
    xki = xk + ai*dx;
    yki = yk + ai*dy;
    fci = feval(fcname,xki);
    fi = fci(1);
    cki = fci(2:end);
    yci = yki - cki;
    ps(i) = fi - tau*sum(log(yki)) + 0.5*bt*(yci'*yci);
end
[psm,ind] = min(ps);
a = aa(ind);
a = min(1,a);