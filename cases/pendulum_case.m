% pendulum_case.m
%
% In this script, we use the data from pendulum simulations. We simulated
% (1) single, (2) double and (3) triple pendulum for chaotic (c) and damped
% (d) cases. These data sets are found in the datasets folder and are
% identified by the name refering to the number of pendulums and case, for
% example a the file '3pd.mat' corresponds to a triple pendulum under a
% damped condition.


clear;  clc;
% Define the number of pendulums (single, double or triple)
pendulums = 2;

% Load the simulation without any noise

data = load([num2str(pendulums), 'p_d.mat']);
y = data.y(:, randperm(length(data.y))); % randomise the order of points.

% Now we rescale and rotate the data set
y = y ./ 5; % rescale the data
R = kron(eye(pendulums), [0, -1; 1, 0]); % rotation matrix
y = R * y; % apply the rotation to the dataset by pi / 2 to center it at 0.

% Add white noise to the data
noise = 0.05; % Set noise level
y = y + noise .* rand(size(y));
N_held = floor(0.25 .* length(y)); % number of points to be held out
y_held = y(:, 1:N_held); % held out data set
y = y(:, N_held + 1:end); % train dataset

% Now we obtain the problem dimensions
M = 2 * pendulums; % Observation dimension
D = pendulums; % Hidden space dimension
N = size(y,2); % Number of data points

% Initialisation step
% All initialisation procedures are functions that can be found in the
% ../init folder. Here are some examples that may be uncommented as
% desired,
% [u_init, A_init, B_init, c_init, lambda2_y_init] = init_random(M, D);
% [u_init, A_init, B_init, c_init, lambda2_y_init] = ...
%     init_true_values(M, D, noise);
scale = 0.25;
[u_init, A_init, B_init, c_init, lambda2_y_init] = ...
    init_perturb(M, D, noise, scale);

% We generate an initial guess for the mean field parameters from u
k1 = repmat(u_init, [1, N]);
k2 = rand(D, N);
m1 = rand(D, N);
m2 = rand(D, N);

% Finally we run the code
[u, A, B, c, lambda2_y, run_stats] = ...
    do_vem(y, u_init, A_init, B_init, c_init, lambda2_y_init,...
    k1, k2, m1, m2, y_held);
       
% Now, we plot the results
% Plotting data (blue) and denoising data (red)
figure
plot(y(1:2:end,:), y(2:2:end,:),'bo')
hold on
plot(y_held(1:2:end,:), y_held(2:2:end,:),'ro')

% Plot key quantities in the model (free energy and prediction error)
figure
subplot(2,1,1)
% obtaining the evolution of the free energy, should go up
fq_hist = run_stats{1}; 
plot(fq_hist);
ylabel('{F}_{q}');
xlabel('Iterations');

subplot(2,1,2)
denoise = run_stats{2};
% obtaining the evolution of the denoising score, should go down
plot(denoise);
ylabel('Prediction RMS');
xlabel('Iterations');

% Plot the color scale
color_scale = [-1.5, 1.5];
figure

% plots for matrix A
subplot(2,2,1)
imagesc(A_init, color_scale)
ylabel('Initial A (M X D)');
subplot(2,2,2)
imagesc(A, color_scale)
ylabel('Final A (M X D)');

% plots for matrix B
subplot(2,2,3)
imagesc(A_init, color_scale)
ylabel('Initial B (M X D)');
subplot(2,2,4)
imagesc(A, color_scale)
ylabel('Final B (M X D)');
