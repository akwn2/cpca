function exitcode = run_shifting(dataset_name, n_latents, noise, save_name, rescale)
%RUN_SHIFTING Runs the case for the given dataset
%   Detailed explanation goes here

if nargin < 2
    n_latents = 1;
end
if nargin < 3
    noise = 1.0;
end
if nargin < 4
    save_name = strcat(dataset_name, '_results_', date());
end
if nargin < 5
    rescale = 0.25;
end

% Fix seed value
% rng(0, 'twister');
randn('seed',0);
rand('seed',0);

% Load dataset
load(dataset_name);
y = y + noise * randn(size(y));

% Rescaling dataset
y = y * rescale;

% Data reordering heur1istic
[y, map] = do_data_reordering(y);

% Separate into held out data and training data
y_held = y(:, 1000:1200);
y = y(:,1:1000);

m_dim = size(y, 1);
n_dim = n_latents;

A = rand(m_dim, n_dim);
B = rand(m_dim, n_dim);

k_prior = rand(n_dim, 1);
m_prior = zeros(n_dim, 1);

u = k_prior .* cos(m_prior);
v = k_prior .* sin(m_prior);

% Initialisatiobn heuristics
lambda2_y = 1/10;
[u, v, A, B, lambda2_y] = init_shifting(y, u, v, A, B, lambda2_y);

save(save_name, 'u', 'v', 'A', 'B', 'lambda2_y', 'y', 'map', 'y_held');

exitcode = 0;
end