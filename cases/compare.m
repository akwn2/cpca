% Run several tests for comparison (each is encapsulated by a for to be
% easier to isolate problems related to indexes and so on).
clear all;
close all;
clc;
% Case code:
% c = chaotic motion data
% d = damped motion data
% m = model-generated data

% rescaling and noise factors
rescaling = 1.;
noise = sqrt(0.1);
savePlots = 0;
useZeroMean = 1;
fixL = 1;

name = ['Annealing init. synthetic data noise ', num2str(noise.^2)];
case_code = ['m']; % case code ['m', 'c', 'd']

for cc = case_code; % case code
    for pp = 1:3 % Pendulum
%         for ll = 1:4 % Latent space size
            ll = pp;
            run_case(strcat(num2str(pp), 'p_', cc,'.mat'), ll, noise, ...
                     strcat(num2str(pp), 'p_', cc,'_out', ...
                        num2str(ll), '.mat'), rescaling, useZeroMean, fixL);
%         end
    end
    make_plots(cc, savePlots, name, useZeroMean);
end

% for cc = case_code; % case code
%     for pp = 1:3 % Pendulum
%         for ll = 1:4 % Latent space size
%             run_shifting(strcat(num2str(pp), 'p_', cc,'.mat'), ll,...
%                         noise, strcat(num2str(pp), 'p_', cc,'_out', ...
%                             num2str(ll), 'shift.mat'), rescaling);
%         end
%     end
%     make_plots_shift(cc, 0);
% end