function [ f, grad ] = optX_dirichlet( M, lambda )
%OPTX_DIRICHLET Calculates dirichlet energy and gradient.

f = @(X) ( lambda * trace(X' * M.L * X) );
grad = @(X) ( lambda * (M.L + M.L') * X );

end

