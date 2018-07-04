addpath('./tools/'),
addpath('/manopt/'), % TODO replace path to manopt if necessary

%% init

M = read_off('./wolf0.off');

params.n_modes = 10; 
params.lambda = 100;
% more parameters available in the functions

%% calculate basis

cmm_normal = cmm(M, params); % this one might be very slow depending on the resolution of the shape / number of modes
cmm_mad = cmm_madmm(M, params);

%% visualize

figure,
for i=1:params.n_modes
    subplot(2,params.n_modes,i), trisurf(M.TRIV, M.VERT(:,1), M.VERT(:,2), M.VERT(:,3), cmm_normal(:,i), 'EdgeAlpha', 0), axis equal, axis off, title('Normal'),
    subplot(2,params.n_modes,i+params.n_modes), trisurf(M.TRIV, M.VERT(:,1), M.VERT(:,2), M.VERT(:,3), cmm_mad(:,i), 'EdgeAlpha', 0), axis equal, axis off, title('MADMM'),
end

colormap(hot),
