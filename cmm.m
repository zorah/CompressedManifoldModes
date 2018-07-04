function [ Phi, Sigma ] = cmm( M, params )
%CMM compute compressed manifold modes of M as described in
% "Compressed Manifold Modes for Mesh Processing" by T. Neumann, 
% K. Varanasi, C. Theobalt, M. magnor and M. Wacker, Eurographics 2014
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
%   params.cmm_tol_abs = 1e-8;
%   params.cmm_tol_rel = 1e-6;
% Outputs:
%   Phi (n x n_modes) - compressed manifold modes
%   Sigma (n_modes) - Dirichlet energy of each function
%
% Written by Zorah Lähner (2016) 
% based on python code by T. Neumann (https://github.com/tneumann/cmm)

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
if ~isfield(params, 'cmm_iterations')
   params.cmm_iterations = 500; 
end
if ~isfield(params, 'cmm_tol_abs')
   params.cmm_tol_abs = 1e-8; 
end
if ~isfield(params, 'cmm_tol_rel')
   params.cmm_tol_rel = 1e-6; 
end

% prepare
[D, L] = calc_cotan_mass(M);
%D = D ./ sum(diag(D));
L = -sparse(L);

% initialize
mu = params.cmm_mu / M.n

Phi = rand(M.n, params.n_modes);
Phi = Phi ./ norm(Phi, 'fro');
E = Phi;
S = Phi;
U_e = zeros(size(Phi));
U_s = zeros(size(Phi));

rho = params.cmm_rho;

D = sqrt(D);
Dinv = diag(1 ./ diag(D));

opts.SYM = true;
opts.POSDEF = true;

% optimize
for i=1:params.cmm_iterations
    i
    % phi step
    Y = 0.5*D*(S - U_s + E - U_e);
    [V, W, U] = svd(Y'*Y);
    Phi = Dinv * (Y*V*(diag(1./sqrt(diag(W))))*U'); % here D actually cancels out... copied from the official code
    
    % E step - TODO use cholesky here to speed up
    A = rho*eye(M.n) - L - L';
    B = rho*(Phi + U_e);
    E_old = E;
    E = linsolve(A, B,opts);
    
    % S step
    v = Phi + U_s;
    max_prox = abs(v) - mu/rho;
    max_prox(max_prox < 0) = 0;
    S_old = S;
    S = sign(v) .* max_prox;
    
    % U step
    U_e = U_e + Phi - E;
    U_s = U_s + Phi - S;
    
    % primal/dual residual & convergence
    if mod(i,params.cmm_check_interval) == 0
        % primal and dual residual
        snorm = sqrt(sum(sum((E - E_old).^2 + (S - S_old).^2)))
        rnorm = sqrt(sum(sum((Phi - E).^2 + (Phi - S).^2)))
        if rnorm > params.cmm_rho_adjust_sensitivity * snorm
            rho = rho * params.cmm_rho_adjust;
            U_s = U_s / params.cmm_rho_adjust;
            U_e = U_e / params.cmm_rho_adjust;
            rho
        elseif snorm > params.cmm_rho_adjust_sensitivity * rnorm
            rho = rho / params.cmm_rho_adjust;
            U_s = U_s * params.cmm_rho_adjust;
            U_e = U_e * params.cmm_rho_adjust;
            rho
        end
        
        % convergence check
        eps_pri = sqrt(size(Phi,2)) * params.cmm_tol_abs + params.cmm_tol_rel * max(max(norm(Phi, 'fro'), norm(E, 'fro')), norm(S, 'fro'));
        eps_dual = sqrt(size(Phi,2)) * params.cmm_tol_abs + params.cmm_tol_rel * norm(rho * U, 'fro');
        if rnorm < eps_pri && snorm < eps_dual
           break; 
        end
    end
end

% calc eigenvalues

Sigma = sum(Phi'*(-L*Phi), 1); % ??

end

