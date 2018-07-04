function [ f, grad ] = optX_reg( rho, Z, U )
%OPTX_REG Calculates regularizer and gradient.

f = @(X) ( rho/2 * sum(sum((X - Z + U).^2)) );
grad = @(X) ( rho * (X - Z + U) );

end

