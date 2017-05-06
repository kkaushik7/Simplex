function [istatus,iB,iN,xB] = simplex_step(A,b,c,iB,iN,xB,irule)
% Take a single simplex method step for the linear program
%
%      min  cx
%      s.t.   Ax=b
%               x>=0,
%
% where A is an (m,n) matrix.
%
% That is, given a basic feasible vector described by the variables iB,iN,xB return the values
% of iB,iN, and xB corresponding to the adjacent basic feasible vector arrived at via a 
% simplex method step.
%
%      Input Parameters:
%
% A - (m,n) constraint matrix
%
% b - (m,1) POSITIVE vector appearing in the constraint equation above
%
% c - (1,n) vector giving the coefficients of the objective function
%
% iB - (1,m) integer vector specifying the indices of the basic
%        variables at the beginning of the simplex step
% iN - (1,n-m) integer vector specifying the indices of the nonbasic 
%         variables at the beginning of the simplex step
% xB - (m,1) vector specifying the values of the basic
%         variables at the beginning of the simplex step
% Binv – (m,m) inverse matrix of the basis B
%
% irule - integer parameter speciying which pivot rule to use:
%       irule = 0 indicates that the smallest coefficient rule should be
%               used
%       irule = 1 indicates that Bland’s rule should be used
%
%             Output Parameters:
%
% istatus - integer parameter reporting on the progress or lake thereof
%               made by this function
%         istatus = 0 indicates normal nondegenerate simplex method step
%               completed
%         istatus = 16 indicates the program is unbounded
%         istatus = -1 indicates an optimal feasible vector has been
%               found
%
% iB - integer vector specifying the m indices of the basic variables
%        after the simplex step 
% iN - integer vector specifying the n-m indices of the nonbasic
%        variables after the simplex step
% xB - vector of length m specifying the values of the basic
%        variables after the simplex step

Binv = inv(A(:,iB));
j = pivotcolumnind(A,c,iB,iN,Binv,irule);
if j == 0 %if no reduced costs are negative
    istatus = -1;
else
    [i,minratio] = pivotrowind(Binv*A(:,j),xB);
    if i == 0 %if no coefficients in column j are positive
        istatus = 16;
    else %if we have a valid pivot entry
        if minratio > 0
            istatus = 0;
        else
            istatus = 0;
        end
        pivotcol = Binv*A(:,j);
        xB = pivot(xB,pivotcol,i); %updates xB
        Binv = pivot(Binv,pivotcol,i); %updates Binv
        iB(i) = j; %update iB
        iN = setdiff(1:size(A,2),iB);
        %update iN
        %for index = 1:length(iN)
        %    if iN(index) == j
        %        iN(index) = i;
        %        break
        %    end
        %end
        
        j = pivotcolumnind(A,c,iB,iN,Binv,irule);
        %if j == 0 %if no reduced costs are negativE
        %    istatus = -1;
        %end
    end
end
end

function [pc] = pivotcolumnind(A,c,iB,iN,Binv,irule)
%Find the appropriate pivot column
%Returns 0 if all reduced costs are nonnegative
pc = 0;
wT = c(iB)*Binv;
if irule == 0
    leastrc = 0;
    for j = iN
        rc = c(j) - wT*A(:,j);
        if rc < leastrc
            pc = j;
            leastrc = rc;
        end
    end
    %OR
    %leastcoeff = 0;
    %for j = iN
    %   if c(j) - wT*A(:,j) < 0 & (c(j)<leastcoeff || pc=0)
    %       pc = j;
    %       leastcoeff = c(j);
    %   end
    %end
elseif irule == 1
    for j = sort(iN)
        if c(j) - wT*A(:,j) < 0
            pc = j;
            break
        end
    end
end
end

function [pr,minratio] = pivotrowind(column,rhs)
%Finds the appropriate pivot row
%Returns 0 if all coefficients in column are nonpositive
pr = 0;
minratio = -1;
    for i = 1:length(column)
        if column(i) > 0
            if rhs(i)/column(i)<minratio || minratio==-1
                pr = i;
                minratio = rhs(i)/column(i);
            end
        end
    end
end

function [newM] = pivot(M,pivotcol,i)
%Returns the matrix obtained by performing on M the elementary matrix
%operation equivalent to pivoting pivotcol at the ith row
newM = M;
m = size(pivotcol);
newM(i,:) = newM(i,:)/pivotcol(i);
for k=1:m
    if k~=i
        newM(k,:)=newM(k,:)-pivotcol(k)*newM(i,:);
    end
end
end