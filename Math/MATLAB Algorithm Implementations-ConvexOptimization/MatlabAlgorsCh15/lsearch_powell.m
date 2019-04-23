% Program: lsearch_powell.m
% Description: Implements Powell's line search method.
% It is used in function sqp_ie_p. 
% Theory: See Sec. 15.3.1.
%======================================================
function a = lsearch_powell(fcname,xk,d_x,muk)
aa = 0:0.01:1;
ps = zeros(101,1);
for i = 1:101,
    ai = aa(i);
    xdi = xk + ai*d_x;
    fci = feval(fcname,xdi);
    fi = fci(1);
    cki = fci(2:end);
    ps(i) = fi - muk'*cki;
end
[psm,ind] = min(ps);
a1 = aa(ind);
ind1 = find(muk <= 1e-5);
s1 = length(ind1);
if s1 == 0,
   a = 0.95*a1;
else
   dk = zeros(s1,1);
   for i = 1:s1,
       for j = 1:101,
           aj = aa(j);
           xdj = xk + aj*d_x;
           fcj = feval(fcname,xdj);
           ckj = fcj(2:end);
           ps(j) = ckj(ind1(i));
       end
       ind2 = find(ps < 0);
       s2 = length(ind2);
       if s2 == 0,
          dk(i) = 1;
       else
          dk(i) = aa(ind2(1)-1);
       end
   end
   a2 = min(dk);
   a = 0.95*min(a1,a2);
end
       
           