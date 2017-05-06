function [L,U] = lu_np(A)
% This function performs LU factorization for
% a matrix A. The function returns the lower
% and upper triangular matrices as separate
% matrices to make checking easier. The matrix
% must be nonsingular and square.
% We will need to know the dimension of the matrix.
dim = size(A,1);
% Initialize L to the identity matrix.
L = eye(dim);
% For each row in A,
for i=1:dim-1,
% For each row under row i,
for j=i+1:dim,
% Check to see if we will encounter divide by zero.
if abs(A(i,i)) <= 1e-14
fprintf('The matrix is singular')
U = NaN;
L = NaN;
return;
end
% Compute the factor.
L(j,i) = A(j,i) / A(i,i);
% Multiply the "nonzero" elements of row i by the
% factor. Subtract this result from the "nonzero"
% elements of row j.
A(j,i+1:dim) = A(j,i+1:dim) - A(i,i+1:dim).*L(j,i);
end
end
% The U factor is the upper triangle of A, with zeros
% in the lower triangle.
U = triu(A);
