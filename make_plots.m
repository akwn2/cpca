function exitcode = make_plots(cc, savePlots, name, useZeroMean, true_A, true_B, true_var)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
close all
if nargin < 2
    savePlots = 0;
end
if nargin < 3
    name = 'unnamed';
end
if nargin < 4
    useZeroMean = 0;
end
% Default true matrices for pendula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    true_A = [1, 0, 0, 0;
            0, 0, 0, 0;
            1, 1, 0, 0;
            0, 0, 0, 0;
            1, 1, 1, 0;
            0, 0, 0, 0;
            1, 1, 1, 1;
            0, 0, 0, 0];

    true_B = [0, 0, 0, 0;
            1, 0, 0, 0;
            0, 0, 0, 0;
            1, 1, 0, 0;
            0, 0, 0, 0;
            1, 1, 1, 0;
            0, 0, 0, 0;
            1, 1, 1, 1];

%     true_A = true_A .* 5;
%     true_B = true_B .* 5;
    
    clim = [-1.5, 1.5];
    k_prior = [ones(1,3), 0.]';
    m_prior = pi / 4 * ones(4,1);

    true_u = k_prior .* cos(m_prior);
    true_v = k_prior .* sin(m_prior);
    true_var = 0.1;
end

% Plot dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:3
    h(ii) = figure;
    load([num2str(ii), 'p_', cc, '_out1.mat']);
    plot(y(1:2:end,:)', y(2:2:end,:)','o');
    title([name,' true latent dimension: ', num2str(ii)]);
    grid on
    legend('first joint', 'second joint', 'third joint');
    if savePlots
%         saveas(gca(), ['./figs/dataset_p_', name, cc, '.fig']);
        hgexport(h(ii), ['./figs/dataset', name,'latent', num2str(ii),'.pdf']);
    end
end


% Plot matrices A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:3
    h(ii) = figure;
    caxis([-5,5])
%         suptitle(strcat('Case: ', num2str(ii), 'p_', cc));
    text(.75, 1.25, strcat(name, ' true latent dimension: ', num2str(ii)));
    ax = zeros(16,1);
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'.mat'));

        ax(4 * jj - 3) = subplot(4, 4, 4 * jj - 3);
        imagesc(A, clim)
        title('Learned A');
        xlabel('D');
        ylabel('M');

        ax(4 * jj - 2) = subplot(4, 4, 4 * jj - 2);
        imagesc(true_A(1:2 * ii, 1:jj), clim)
        title('True A');
        xlabel('D');
        ylabel('M');

        ax(4 * jj - 1) = subplot(4, 4, 4 * jj - 1);
        imagesc(B, clim)
        title('Learned B');
        xlabel('D');
        ylabel('M');

        ax(4 * jj) = subplot(4, 4, 4 * jj);
        imagesc(true_B(1:2 * ii, 1:jj), clim)
        title('True B');
        xlabel('D');
        ylabel('M');
    end
    axis tight
    gcf();
    cb = colorbar;
    set(cb, 'Position', [.90 .125 .0525 .8])
    for jj=1:16
          pos=get(ax(jj), 'Position');
          set(ax(jj), 'Position', [pos(1) pos(2) 0.75 * pos(3) pos(4)]);
    end
    
    if savePlots
%         saveas(gca(), ['./figs/matrices_case', num2str(ii), 'p_', name, cc, '.fig']);
        hgexport(h(ii), ['./figs/matrices', name,' latent ', num2str(ii), '.pdf']);
    end
end



% Plot variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
score = ones(4,1);
for ii = 1:3
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'.mat'));
        score(jj) = 1 / lambda2_y;
    end
    subplot(3, 1, ii)
    plot(score, '-*b');
    hold on
    plot(true_var .* ones(4,1),'-+r');
    legend('Learned','True')
    title(strcat(name,' (x,y) pair: ', num2str(ii)));
    ylim([0, max([true_var, score']) + 1.]);
    ylabel('Variance');
    xlabel('Latent space size');
end
if savePlots
%     saveas(gca(), ['./figs/var_p_', name, cc, '.fig']);
    hgexport(h, ['./figs/var', name,'.pdf']);
end


% Plot prior concentration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
learned_u = zeros(3,4);
for ii = 1:3
    for jj = 1:4
        load([num2str(ii),'p_',cc,'_out',num2str(jj),'.mat']);
        learned_u(ii,1:jj) = u(1:jj);
        subplot(3, 4, 4 * (ii - 1) + jj)
        plot(learned_u(ii,:), '-*');
        hold on
        plot(true_u,'-+r');
        xlabel(strcat('N = ', num2str(jj)));
        title(strcat(name,'', num2str(ii), ' pendulum ', cc));
        legend('Learned', 'True')
        ylabel('Prior concentration');
    end
end
if savePlots
%     saveas(gca(), ['./figs/concentration', name, cc, '.fig']);
    hgexport(h, ['./figs/concentration', name,'.pdf']);
end


% Plot error score (denoising)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
rms = zeros(3,4);
true_rms = zeros(3,4);
for ii = 1:3
    for jj = 1:4
        load([num2str(ii),'p_',cc,'_out',num2str(jj),'.mat']);
        rms(ii, jj) = run_stats{3}(end);
        
        % Get ground truth rms
        if useZeroMean
            [m_dim, n_dim] = size(A);
            d_held = size(y_held, 2);

            model_pars = zm_pack_model_pars(true_u(1:n_dim), ...
                true_A(1:m_dim, 1:n_dim), true_B(1:m_dim, 1:n_dim), 1/true_var);

            pars = pack_pars(y_held, rand(n_dim, d_held), rand(n_dim, d_held), ...
                rand(n_dim, d_held), rand(n_dim, d_held), rand(n_dim, d_held), ...
                n_dim, m_dim, d_held);

            true_rms(ii, jj) = zm_get_validation_error(model_pars, pars, y_held, d_held);
        else
            [m_dim, n_dim] = size(A);
            d_held = size(y_held, 2);

            model_pars = pack_model_pars(true_u(size(u)), true_v(size(v)), ...
                true_A(1:m_dim, 1:n_dim), true_B(1:m_dim, 1:n_dim), 1/true_var);

            pars = pack_pars(y_held, rand(n_dim, d_held), rand(n_dim, d_held), ...
                rand(n_dim, d_held), rand(n_dim, d_held), rand(n_dim, d_held), ...
                n_dim, m_dim, d_held);

            true_rms(ii, jj) = get_validation_error(model_pars, pars, y_held, d_held);
        end
    end
    subplot(3,1,ii)
    plot(rms(ii,:), '-*b')
    hold on
    plot(true_rms(ii,:), '-+r')
    title([name, ' true latent dimension: ', num2str(ii)]);
    legend('Learned',' True')
    xlabel('Latent space dimension');
    ylabel('RMS');
end
if savePlots
%     saveas(gca(), ['./figs/rms', name, cc, '.fig']);
    hgexport(h, ['./figs/rms', name,'.pdf']);
end


exitcode = 0;
end