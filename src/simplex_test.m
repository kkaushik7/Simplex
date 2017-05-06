function [X,eta,isfeasible,isoptimal,zN]=simplex_test(A,b,c,iB,xB)
% function [X,eta,isfeasible,isoptimal,zN]=simplex_test(A,b,c,iB,xB)
% 
% Test the feasibility and optimality of a basic vector for the
% linear program:
%
%     min:          cx
%     subject to:   Ax=b                                    (1)
%                   x>= 0
%
% where A is an (m,n) matrix of rank m.
%
%                     Input Parameters:
%
% A - (m,n) constraint matrix
% b - (m,1) vector appearing in the constraint equation above
% c - (1,n) vector giving the coefficients of the objective function
%
% iB - (1,m) integer vector specifying the indices of the basic
%      variables at the beginning of the simplex step
% xB - (m,1) vector specifying the values of the basic 
%      variables at the beginning of the simplex step
%
%                     Output Parameters:
%
%  X - vector of length n which contains both the basic and
%    nonbasic values
%  eta - value of the objective function at X
%  isfeasible - integer parameter indicating whether the given basic
%    vector is feasible or not;
%
%      isfeasible = 0  means the vector is infeasible
%      isfeasible = 1  means the vector is feasible
%
%  zN - the value of the dual basic variables corresponding to the
%    given basic feasible vector
%

[m,n] = size(A);
eps = 1.0e-12;

% initialize the output variables
X=[];
eta=[];
isfeasible=[];
isoptimal=[];

% form the entire vector X
X = zeros(n,1);
X(iB)=xB;
eta=c*X;

% check for feasiblity
err=norm(A*X-b);

isfeasible=0;
if (err < eps) && isempty(find(X<-eps))
   isfeasible=1;
end

% form the integer vector iN
iN=1:n;
iN(iB)=0;
iN=iN(find(iN));

% check for optimality
c1=c(iB);
c2=c(iN);

B = A(:,iB);
N = A(:,iN);

Binv = inv(B);
ctilde = c2-c1*Binv*N;

if ~isempty(find(ctilde<0)) 
   isoptimal=0;
else
   isoptimal=1;
end

zN=ctilde;
