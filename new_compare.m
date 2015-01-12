% Run several tests for comparison (each is encapsulated by a for to be
% easier to isolate problems related to indexes and so on).
clear all;
close all;
clc;
% Case code:
% c = chaotic motion data
% d = damped motion data
% m = model-generated data

%initialisation types:
% 'true' - using true parameters
% 'random' - randomly chosen parameters
% 'annealing' - via annealing scheme
% 'bootstrap' - by bootstrapping

% rescaling and noise factors
rescale = 1.;
sigma2 = 0.025;
noise = sqrt(sigma2);
savePlots = 1;
useZeroMean = 1;
reorderData = 0;

case_code = ['m']; % ['m', 'c', 'd']
initType = {'true','random','annealing'}; %{'true','random','annealing','bootstrap'};

for cc = case_code; % case code
    for pp = 1:3
        for ii = initType
            name = ['Case=', cc,', '...
                    'Init=', strjoin(ii),', '...
                    'Variance=', num2str(noise.^2)];

            datasetName = [num2str(pp), 'p_', cc,'.mat'];
            
            saveName = [num2str(pp), 'p_', cc, '_', strjoin(ii), '_out.mat'];

            new_run_case(datasetName, saveName, useZeroMean, ...
                reorderData, strjoin(ii), pp, noise, rescale);
            
            new_make_plots(cc, pp, strjoin(ii), noise, useZeroMean, ...
                savePlots);
        end
    end
end