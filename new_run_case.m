function exitcode = new_run_case(datasetName, saveName, useZeroMean, ...
    reorderData, initType, n_latents, noise, rescale)
%NEW_RUN_CASE Runs the case for the given dataset
%   Detailed explanation goes here

if nargin < 2
    saveResults = 0;
else
    saveResults = 1;
end

if nargin < 3
    useZeroMean = 0;
end

if nargin < 4
    reorderData = 0;
end

if nargin < 5
    initType = '';
end

if nargin < 6
    n_latents = 1;
end

if nargin < 7
    noise = 1.0;
end

if nargin < 8
    rescale = 0.25;
end


% Fix seed value
try
    rng(0, 'twister');
catch
    randn('seed',0);
    rand('seed',0);
end


% Load, add noise and rescale dataset. Then separate into held out data
% and training data
load(datasetName);
y = y + noise * randn(size(y));
y = y * rescale;
y = y(:,randperm(1200));
y_held = y(:, 1000:1200);
y = y(:,1:1000);

[m_dim, d_pts] = size(y);
n_dim = n_latents;

% Data reordering heur1istic
if reorderData
    [y, map] = do_data_reordering(y);
else
    map = 1:m_dim;
end

% Initialisation schemes
switch initType
    case 'annealing'
        [u, v, A, B, ~] = init_randomly(m_dim, n_dim, useZeroMean);
        [u, v, A, B, lambda2_y] = init_annealing(y, u, v, A, B, ...
            y_held, useZeroMean);
    case 'bootstrap'
        [u, v, A, B, lambda2_y] = init_randomly(m_dim, n_dim, useZeroMean);
        [u, v, A, B, lambda2_y] = init_bootstrap(y, u, v, A, B, ...
            lambda2_y, y_held, useZeroMean);
    case 'true'
        [u, v, A, B, lambda2_y] = init_true_val(m_dim, n_dim, ...
            useZeroMean, noise);
    case 'random'
        [u, v, A, B, lambda2_y] = init_randomly(m_dim, n_dim, useZeroMean);
    otherwise
        fprintf('unrecognized case, aborting!');
        exitcode = 1;
        return
end
fprintf('\n\nInitial guesses\n\n');
disp(u)
disp(v)
disp(A)
disp(B)
disp(lambda2_y)

A = A * rescale;
B = B * rescale;

fprintf('\n\nMain iteration\n\n');

k1 = repmat(abs(u + 1.i * v), [1, d_pts]);
k2 = small_rand(n_dim, d_pts);
m1 = repmat(angle(u + 1.i * v), [1, d_pts]);
m2 = small_rand(n_dim, d_pts);

if useZeroMean
    [u, A, B, lambda2_y, run_stats] = ...
        zm_do_vem(y, u, A, B, lambda2_y, k1, k2, m1, m2, y_held);
else
    [u, v, A, B, lambda2_y, run_stats] = ...
        do_vem(y, u, v, A, B, lambda2_y, k1, k2, m1, m2, y_held);
end

A = A ./ rescale;
B = B ./ rescale;
y = y ./ rescale;
y_held = y_held ./ rescale;

if saveResults
    save(saveName, 'u', 'v', 'A', 'B', 'lambda2_y', 'y', 'run_stats', 'map', 'y_held');
end


exitcode = 0;
end