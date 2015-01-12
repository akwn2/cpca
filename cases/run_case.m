function failed = run_case(dataset_name, n_latents, noise, save_name, rescale, useZeroMean, fixL)
%RUN_CASE Runs the case for the given dataset
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
if nargin < 6
    useZeroMean = 0;
end
if nargin < 7
    fixL = 0;
end

% Fix seed value
%rng(0, 'twister');
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
y = y(:,randperm(1200));
y_held = y(:, 1000:1200);
y = y(:,1:1000);

m_dim = size(y, 1);
n_dim = n_latents;

weight = 1;
 
A = rand(m_dim, n_dim) / weight;
B = rand(m_dim, n_dim) / weight;

k_prior = rand(n_dim, 1) / weight;
m_prior = zeros(n_dim, 1) / weight;

if useZeroMean
    u = k_prior;
else
    u = k_prior .* cos(m_prior);
    v = k_prior .* sin(m_prior);
end
% Initialisatiobn heuristics
sigma2 = 100;
lambda2_y = 1/sigma2;
% [u, v, A, B, lambda2_y] = init_shifting(y, u, v, A, B, lambda2_y);
% [u, v, A, B, lambda2_y] = init_bootstrap(y, u, v, A, B, lambda2_y);
% [u, v, A, B, lambda2_y] = init_sequential(y, u, v, A, B, lambda2_y);

fprintf('\n\nMain iteration\n\n');

d_pts = size(y, 2);

k1 = repmat(k_prior, [1, d_pts]);
k2 = rand(n_dim, d_pts);
m1 = repmat(m_prior, [1, d_pts]);
m2 = rand(n_dim, d_pts);

if useZeroMean
    if fixL   
        for ll = 1:4
            lambda2_y = lambda2_y / ll;
            [u, A, B, lambda2_y, run_stats] = ...
                fixL_do_vem(y, u, A, B, lambda2_y, k1, k2, m1, m2, y_held);
            v = zeros(size(u));
        end
    end
    [u, A, B, lambda2_y, run_stats] = ...
        zm_do_vem(y, u, A, B, lambda2_y, k1, k2, m1, m2, y_held);
    v = zeros(size(u));
else
    [u, v, A, B, lambda2_y, run_stats] = ...
        do_vem(y, u, v, A, B, lambda2_y, k1, k2, m1, m2, y_held);
end
%     [y_model, y_model_noiseless] = generate_data_from_model(u, A, B, lambda2_y, d_pts);

A = 1/rescale * A;
B = 1/rescale * B;
y = 1/rescale * y;
y_held = 1/rescale * y_held;
%     y_model = 1/rescale * y_model;
%     y_model_noiseless = 1/rescale * y_model_noiseless;

save(save_name, 'u', 'v', 'A', 'B', 'lambda2_y', 'y', 'run_stats', 'map', 'y_held');

failed = 0;
% catch
%     failed = 1;
% end
end
