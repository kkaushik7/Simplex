function [istatus,iB,iN,xB] = simplex_init(A,b,c)
%
% Attempt to find a basic feasible vector for the linear program
%
%     min    cx
%     s.t.     Ax=b
%               x>=0,
%
% where A is a (m,n) matrix.
%
%         Input Parameters:
%
% A - (m,n) constraint matrix
% b - (m,1) vector appearing in the constraint equation above
% c - (1,n) vector giving the coefficients of the objective function
%
%        Output Parameters:
%
% istatus - integer parameter reporting the result of the initialization procedure
%       istatus = 0 indicates a basic feasible vector was found
%       istatus = 4 indicates that the initialization procedure failed
%       istatus = 16 indicates that the problem is infeasible
%
% iB - integer vector of length m specifying the indices of the basic variables
% iN - integer vector of length n-m specying the indices of the nonbasic variables

% xB - vector of length m specifying the values of the basic variables
%

%Find elements of b which are negative
[m,n] = size(A);
E = eye(m);
for index = 1:m
    if b(index) < 0
        E(index,index) = -1;
    end
end

A1 = [E*A eye(m)];
b1 = E*b;
c1 = [zeros(1,n) ones(1,m)];
iB1 = n+1:n+m;
iN1 = 1:n;
Binv1 = eye(m);
xB1 = b1;
irule=1;

iB = [];
iN = [];
xB = [];
istatus = 0;
while istatus ~= -1
    [istatus,iB1,iN1,xB1]=simplex_step(A1,b1,c1,iB1,iN1,xB1,irule);
end
if c1(iB1)*xB1 == 0
    if max(iB1) > n
        istatus = 4;
    else
        istatus = 0;
        iB = iB1;
        iN = setdiff(1:n,iB);
        xB = xB1;
        
    end
    
else
    istatus = 16;
    
  
end