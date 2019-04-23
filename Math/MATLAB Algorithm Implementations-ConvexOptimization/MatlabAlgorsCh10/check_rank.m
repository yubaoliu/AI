% Program: check_rank.m
% Description: Checks the equality constraints for a linear system 
% Ax = b and evaluates the ranks of matrices A and 
% [A b].
% Theory: See Practical Optimization Sec. 10.2.2.
% Input:
% A: coefficient matrix of linear system Ax = b
% b: vector b in linear system Ax = b
% Output:
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
% [p,p1,p2] = check_rank(A,b)
% ==========================================================
function [p,p1,p2] = check_rank(A,b)
 disp(' ')
 disp('Program check_rank.m')
sz = size(A);
p = sz(1);
p1 = rank(A);
p2 = rank([A b]);
if p1 ~= p2,
disp('The linear equality constraints are inconsistent.')
elseif p1 < p,
disp('The linear equality constraints are consistent but redundant.')
else
disp('The linear equality constraints are consistent without redundancy.')
end