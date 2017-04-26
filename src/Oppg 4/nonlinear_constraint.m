function [ c, ceq ] = nonlinear_constraint( z ) 
%   Implements the nonlinear equality constraint c(z_k) <= e_k
% x1 = lambda_k, x5 = e_k
global N mx beta alpha lambda_t 
lambda_k = z(1:mx:N*mx);
e_k = z(5:mx:N*mx);

c = alpha*exp(-beta*(lambda_k - lambda_t).^2) - e_k;
ceq = [];
end

