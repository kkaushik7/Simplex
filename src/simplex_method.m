function [istatus,X,eta,iB,iN,xB] = simplex_method(A,b,c,irule)
%
% Find a basic optimal solution for the linear program
%
%    min   cx
%    s.t.    Ax=b
%               x>=0,
%
% where A is an (m,n) matrix.
%
%         Input Parameters:
%
% A - (m,n) constraint matrix
% b - (m,1) a vector appearing in the constraint equation above
% c - (1,n) vector giving the coefficients of the objective function
%
% irule - integer parameter speciying which pivot rule to use:
%      irule = 0 indicates that the smallest coefficient rule should be used
%      irule = 1 indicates that Bland’s rule should be used
%
%        Output Parameters:
%

% istatus - integer parameter reporting the results obtained by this function
%       istatus = 0 indicates normal completion (i.e., a solution has been found and reported)
%       istatus = 4 indicates the program is infeasible
%       istatus = 16 indicates the program is feasible but our initialization procedure has failed 
%       istatus = 32 indicates that the program is unbounded
%
% X - vector of length n specifying the solution 
% eta - the minimum value of the objective function
% iB - integer vector specifying the m indices of the basic variables after the simplex step
% iN - integer vector specifying the n-m indices of the nonbasic variables after the simplex step
% xB - vector of length m specifying the values of the basic variables after the simplex step
%
[istatus,iB,iN,xB] = simplex_init(A,b,c);

if istatus == 4 || istatus == 16
    X = NaN;
    eta = NaN;
    if istatus == 4
        istatus = 16;
    else
        istatus = 4;
    end
else
    [m,n] = size(A);
    Binv = A(:,iB)\eye(m);
    while istatus ==0
        [istatus,iB,iN,xB] = simplex_step(A,b,c,iB,iN,xB,irule);
    end
    if istatus == -1
        X = zeros(n,1);
        X(iB) = xB;
        eta = c(iB)*xB;
        istatus = 0;
    else
        X = NaN;
        eta = NaN;
        istatus = 32;
    end
end
end