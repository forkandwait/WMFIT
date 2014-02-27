function [S, iterations, errors]  = ipf(T, R, C, epsilon, max_it)
%
% file:      	ipf.m, (c) Matthew Roughan, Tue Feb 17 2009
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
% Iterative Proportional Fitting (IPF) 
%     There are a couple of ways to use IPF, but in this case the aim is to find a matrix close
%     to the input T that satisfies contraints on its row sums and column sums (its marginals). In
%     particular we aim to find S with row sums equal to R, and column sums equal to C.
%
% Note that because sum(R) and sum(T) are both sum(sum(T)), we must have sum(R)=sum(C) in the
% input parameters. Also, the row and col sums R and C must be >= 0, or we end up dividing by
% zero. Also, T must not have any zero row or column sums.
%
% The algorithm essentially works by going through the constraints (each row, and then each
% col constraint) in turn, and scaling the appropriate row or col of the matrix to ensure
% that it satisfies the constraint. The magic of IPF is that given appropriate inputs, it
% will always converge to a solution that satisfies the constraints.
%
% INPUTS:
%     T       : nxm non-negative input matrix
%     R       : nx1 vector of required row sums
%     C       : 1xm vector of required column sums
%     epsilon : optional precision for convergence (default 1e-3)
%     max_it  : optional maximum number of iterations (default 100)
%
% OUTPUTS:
%     S:          IPF'd version of T
%     iterations: number of iterations until convergence
%     errors:     (1 x iterations) vector of errors, so that we can examine how quickly it
%                 converged 
% 
% see ipf_test.m for some examples of how to call the function
%
% References: easiest place to look is Wikipedia
%      http://en.wikipedia.org/wiki/Iterative_proportional_fitting
% where the algorithm implemented here corresponds to algorithm 3 (RAS)
%     

% test inputs
st = size(T);
sr = size(R);
sc = size(C);
if (st(1) ~= sr(1))
  error('T or R is the wrong size.');
end
if (st(2) ~= sc(2))
  error('T or C is the wrong size.');
end
if (sum(sum(T < 0)) > 0)
  error('T must not contain negative values.');
end
if (sum(R<=0)>0 || sum(C<=0)>0)
  error('R and C must be positive.');
end
if (abs(sum(R)-sum(C)) > 1e-6*sum(R))
  error('R and C must have the same sum.');
end
if (sum(sum(T,1)==0)>0 ||  sum(sum(T,2)==0)>0)
  error('row and column sums of T must be non-zero.');
end
if (nargin < 4)
  epsilon = 1e-3;
end
if (nargin < 5)
  max_it = 100;
end
n = st(1);
m = st(2);

% iterations
total_error = sum(abs(sum(T,2)-R)) + sum(abs(sum(T,1)-C));
S = T;
count = 1;
errors(count) = total_error;

while (total_error > epsilon && count < max_it)
	% scale rows 
	%    could probably write smarter matlab here to vectorize the whole operation using diag's
	%    but the number of operations would go up, so it might not be faster anyway
	for i=1:n
		S(i,:) = S(i,:) * R(i) / sum(S(i,:));
	end
	
	% scale columns
	for j=1:m
		S(:,j) = S(:,j) * C(j) / sum(S(:,j));
	end
	
	% test errors
	total_error = sum(abs(sum(S,2)-R)) + sum(abs(sum(S,1)-C));
	count = count + 1;
	errors(count) = total_error;
end % while

iterations = count;
