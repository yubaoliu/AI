% Program: check_rank_m.m
% Description: Checks the linear equality constraints for 
% a linear system Ax = b and evaluates the ranks of matrices A
% and [A b]. In the case where the equations in the system are 
% consistent but redundant, the function produces an equivalent 
% but reduced-size linear system in which the redundancy 
% is removed.
% Theory: See Practical Optimization Sec. 10.2.2.
% Input:
% A: coefficient matrix of linear system Ax = b
% b: vector b in linear system Ax = b
% Output:
% Ah: (a) Ah is the same as A if system Ax = b is 
%         consistent with no redundancy;
%     (b) Ah is the coefficient matrix of the reduced system 
%         Ah*x = bh that is equivalent to system Ax = b
%         without redundancy, if Ax = b is a consistent 
%         linear system with redundancy.
%     (c) Ah is empty if Ax = b is inconsistent.
% bh: (a) bh is the same as b if system Ax = b is 
%         consistent and no redundancy;
%     (b) bh is the constantvector in reduced system 
%         Ah*x = bh that is equivalent to system Ax = b
%         wothout redundancy, if Ax = b is a consistent 
%         linear system with redundancy.
%     (c) bh is empty if Ax = b is inconsistent.
% p: number of equations in Ax = b
% p1: rank of matrix A
% p2: rank of matrix [A b]
% Example:
% Find the ranks of A and [A b] in the linear system
% 2*x1 - 5*x2 = -6
% -x1 + 2.5*x2 = 3
% Solution:
% Execute the command
% A = [2 -5; -1 2.5]
% b = [-6; 3]
% [Ah,bh,p,p1,p2] = check_rank_m(A,b)
function [Ah,bh,p,p1,p2] = check_rank_m(A,b)
 disp(' ')
 disp('Program check_rank.m')
sz = size(A);
p = sz(1);
p1 = rank(A);
p2 = rank([A b]);
if p1 ~= p2,
   disp('The linear equality constraints are inconsistent.')
   Ah = [];
   bh = [];
elseif p1 < p,
   disp('The linear equality constraints are consistent but redundant.')
   [U,S,V] = svd(A);
   S1 = S(1:p1,1:p1);
   Ah = S1*(V(:,1:p1))'; 
   b1 = U'*b;
   bh = b1(1:p1);
else
   disp('The linear equality constraints are consistent without redundancy.')
   Ah = A;
   bh = b;
end