function [ Phi, Sigma ] = cmm_madmm( M, params )
%CMM_MADMM compute compressed manifold modes of M as described in
% "Compressed Manifold Modes for Mesh Processing" by T. Neumann, 
% K. Varanasi, C. Theobalt, M. magnor and M. Wacker, Eurographics 2014 [1]
% and "MADMM: a generic algorithm for non-smooth optimization on
% manifolds" by A. Kovnatsky, M. M. Bronstein and K. Glashoff, ECCV 2016
% [2]
% The code combines the optimization of the Stiefel manifold [2] with
% primal/dual residuals and convergence check from [1]. 
% Inputs:
%   M - struct with M.VERT (n x 3), M.TRIV (m x 3), M.n (number of vertices)
%  Possible params fields and default values:
%   params.n_modes = 10;
%   params.cmm_mu = 1;
%   params.cmm_rho = 1;
%   params.cmm_check_interval = 10;
%   params.cmm_rho_adjust_sensitivity = 5;
%   params.cmm_rho_adjust = 2;
%   params.cmm_iterations = 500;
%   params.cmm_lambda = 1;
%   params.cmm_tol_abs = 1e-8;
%   params.cmm_tol_rel = 1e-6;
%   params.silent = false;
%   params.manopt_maxiter = 3;
%   params.manopt_tolgradnorm = 1e-6;
%   params.manopt_manstepsize = 1e-6;
%   params.manopt_verbosity = 0;
% Outputs:
%   Phi (n x n_modes) - compressed manifold modes
%   Sigma (n_modes) - Dirichlet energy of each function
%
% Written by Zorah Lähner (2016)

% set parameters
if ~isfield(params, 'n_modes')
   params.n_modes = 10; 
end
if ~isfield(params, 'cmm_mu')
   params.cmm_mu = 1; 
end
if ~isfield(params, 'cmm_rho')
   params.cmm_rho = 1; 
end
if ~isfield(params, 'cmm_check_interval')
   params.cmm_check_interval = 10; 
end
if ~isfield(params, 'cmm_rho_adjust_sensitivity')
   params.cmm_rho_adjust_sensitivity = 5; 
end
if ~isfield(params, 'cmm_rho_adjust')
   params.cmm_rho_adjust = 2; 
end
if ~isfield(params, 'cmm_lambda')
   params.cmm_lambda = 1; 
end
if ~isfield(params, 'cmm_iterations')
   params.cmm_iterations = 500; 
end
if ~isfield(params, 'cmm_tol_abs')
   params.cmm_tol_abs = 1e-8; 
end
if ~isfield(params, 'cmm_tol_rel')
   params.cmm_tol_rel = 1e-6; 
end
if ~isfield(params, 'silent')
    params.silent = false;
end
% set manopt options
if ~isfield(params, 'manopt_maxiter')
    options.maxiter = 3;
else
    options.maxiter = params.manopt_maxiter;
end
if ~isfield(params, 'manopt_tolgradnorm')
    options.tolgradnorm = 1e-6;
else
    options.tolgradnorm = params.manopt_tolgradnorm;
end
if ~isfield(params, 'manopt_manstepsize')
    options.minstepsize = 1e-6;
else
    options.minstepsize = params.manopt_minstepsize;
end
if ~isfield(params, 'manopt_verbosity')
    options.verbosity = 0;
else
    options.verbosity = params.manopt_verbosity;
end

options.stopfun = @mystopfun;
    function stopnow = mystopfun(problem, x, info, last)
        stopnow = (last >= 3 && (info(last-2).cost - info(last).cost)/info(last).cost < 1e-8);
    end

% prepare
[M.A, M.S] = calc_cotan_mass(M);
M.S = M.S;
M.L = sparse(M.S);

mu = params.cmm_mu  / M.n;
rho = params.cmm_rho;

% stiefel factory for X step
X_manifold = stiefelfactory(M.n,params.n_modes);
X_problem = {};
X_problem.M = X_manifold;

[data_f, data_grad] = optX_dirichlet(M, params.cmm_lambda);

%[evecs, evals] = eigs(M.L, 10);

% init
%X = orth(rand(M.n, params.n_modes));
% Z = orth(rand(M.n, params.n_modes));
[X,~] = svd(randn(M.n,params.n_modes),0);
Z = X;
U = X - Z;

for k=1:params.cmm_iterations
    if ~params.silent 
        k 
    end
    % X step
    manifold = stiefelfactory(M.n,params.n_modes);
    problem = {};
    problem.M = manifold;
    
    [reg_f, reg_grad] = optX_reg(rho, Z, U);
    problem.cost = @(X)  (data_f(X) + reg_f(X) );
    problem.egrad = @(X) ( data_grad(X) + reg_grad(X) );
    
    %checkgradient(problem);
    
    %[X, cost, ~, ~] = trustregions(problem, X, options);
    [X, cost, ~, ~] = conjugategradient(problem, X, options);
    if ~params.silent 
        cost
    end
    
    % Z step
    v = X + U;
    Z_old = Z;
    Z = sign(v) .* max(0, abs(v) - mu/rho );
    
    % U update
    U = U + X - Z;
    
    % primal/dual residual & convergence
    if mod(k,params.cmm_check_interval) == 0
        % primal and dual residual
        snorm = sqrt(sum(sum((Z - Z_old).^2)));
        rnorm = sqrt(sum(sum((X - Z).^2)));
        if rnorm > params.cmm_rho_adjust_sensitivity * snorm
            rho = rho * params.cmm_rho_adjust;
            U = U / params.cmm_rho_adjust;
            if ~params.silent
                rho
            end
        elseif snorm > params.cmm_rho_adjust_sensitivity * rnorm
            rho = rho / params.cmm_rho_adjust;
            U = U * params.cmm_rho_adjust;
            if ~params.silent
                rho
            end
        end
        
        % convergence check
        eps_pri = sqrt(size(X,2)) * params.cmm_tol_abs + params.cmm_tol_rel * max(max(norm(X, 'fro'), norm(Z, 'fro')));
        eps_dual = sqrt(size(X,2)) * params.cmm_tol_abs + params.cmm_tol_rel * norm(rho * U, 'fro');
        if rnorm < eps_pri && snorm < eps_dual
           break; 
        end
    end
    
end

Phi = X;
Sigma = sum(Phi'*(M.S*Phi), 1); 

end

