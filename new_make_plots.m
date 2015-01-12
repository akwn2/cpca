function exitcode = new_make_plots(cc, jj, init, noise, useZeroMean, savePlots, clim)
%NEW_MAKE_PLOTS Makes plots for cases ran
%   Detailed explanation goes here
close all

if nargin < 5
    useZeroMean = 0;
end
if nargin < 6
    savePlots = 0;
end
if nargin < 7
    clim = [-1.5, 1.5];
end
load([num2str(jj), 'p_', cc, '_', init, '_out.mat']);
[m_dim, n_dim] = size(A);
d_pts = size(y, 2);
d_held = size(y_held, 2);

[true_u, true_v, true_A, true_B, true_lambda2_y] = ...
    init_true_val(m_dim, n_dim, useZeroMean, noise);

% Plot dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
plot(y(1:2:end,:)', y(2:2:end,:)','o');

name = ['Case=', cc,', '...
        'Init=', init,', '...
        'Variance=', num2str(noise.^2),', '...
        'N=', num2str(d_pts), ', '...
        'M=', num2str(n_dim), ', '...
        'D=', num2str(m_dim), ];

title(name);
grid on
legend('first joint', 'second joint', 'third joint');

if savePlots
    saveas(gca(), ['dataset', name, num2str(jj), '.pdf']);
%     hgexport(h, ['dataset', name, num2str(jj), '.pdf']);
end


% Plot matrices A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
caxis([-5,5])
try
    suptitle(name);
catch
    text(.75, 1.25, name);
end

ax(1) = subplot(4, 1, 1);
imagesc(A, clim)
title('Learned A');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

ax(2) = subplot(4, 1, 2);
imagesc(true_A, clim)
title('True A');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

ax(3) = subplot(4, 1, 3);
imagesc(B, clim)
title('Learned B');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

ax(4) = subplot(4, 1, 4);
imagesc(true_B, clim)
title('True B');
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

axis tight
gcf();
cb = colorbar;
set(cb, 'Position', [.90 .125 .0525 .8])

for mm=1:4
      pos=get(ax(mm), 'Position');
      set(ax(mm), 'Position', 0.95 * [pos(1) pos(2) 0.9 * pos(3) pos(4)]);
end

if savePlots
    saveas(gca(), ['matrices', name, num2str(jj), '.pdf']);
%     hgexport(h, ['matrices', name, num2str(jj), '.pdf']);
end


% % 
% % % Plot variances
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = figure;
% % score = ones(4,1);
% % for jj = 1:3
% %     for jj = 1:4
% %         load(strcat(num2str(jj),'p_',cc,'_out',num2str(jj),'.mat'));
% %         score(jj) = 1 / lambda2_y;
% %     end
% %     subplot(3, 1, jj)
% %     plot(score, '-*b');
% %     hold on
% %     plot(true_var .* ones(4,1),'-+r');
% %     legend('Learned','True')
% %     title(strcat(name,' (x,y) pair: ', num2str(jj)));
% %     ylim([0, max([true_var, score']) + 1.]);
% %     ylabel('Variance');
% %     xlabel('Latent space size');
% % end
% % if savePlots
% % %     saveas(gca(), ['var_p_', name, cc, '.fig']);
% %     hgexport(h, ['var', name,'.pdf']);
% % end
% % 
% % 
% % % Plot prior concentration
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % h = figure;
% % learned_u = zeros(3,4);
% % for jj = 1:3
% %     for jj = 1:4
% %         load([num2str(jj),'p_',cc,'_out',num2str(jj),'.mat']);
% %         learned_u(jj,1:jj) = u(1:jj);
% %         subplot(3, 4, 4 * (jj - 1) + jj)
% %         plot(learned_u(jj,:), '-*');
% %         hold on
% %         plot(true_u,'-+r');
% %         xlabel(strcat('N = ', num2str(jj)));
% %         title(strcat(name,'', num2str(jj), ' pendulum ', cc));
% %         legend('Learned', 'True')
% %         ylabel('Prior concentration');
% %     end
% % end
% % if savePlots
% % %     saveas(gca(), ['concentration', name, cc, '.fig']);
% %     hgexport(h, ['concentration', name,'.pdf']);
% % end


% Plot error score (denoising)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
rms = zeros(3,4);
true_rms = zeros(3,4);
        
if useZeroMean
    pars = pack_pars(y_held, small_rand(n_dim, d_held), ...
        small_rand(n_dim, d_held), small_rand(n_dim, d_held), ...
        small_rand(n_dim, d_held), small_rand(n_dim, d_held), ...
        n_dim, m_dim, d_held);

    % Get ground paramters rms
    model_pars = zm_pack_model_pars(u, A, B, true_lambda2_y);

    rms = zm_get_validation_error(model_pars, pars, ...
        y_held, d_held);

    % Get ground truth rms
    model_pars = zm_pack_model_pars(true_u, ...
        true_A, true_B, true_lambda2_y);

    true_rms = zm_get_validation_error(model_pars, pars, ...
        y_held, d_held);
else
    pars = pack_pars(y_held, small_rand(n_dim, d_held), ...
        small_rand(n_dim, d_held), small_rand(n_dim, d_held), ...
        small_rand(n_dim, d_held), small_rand(n_dim, d_held), ...
        n_dim, m_dim, d_held);

    % Get ground paramters rms
    model_pars = pack_model_pars(u, v, A, B, true_lambda2_y);

    rms = get_validation_error(model_pars, pars, y_held, d_held);

    % Get ground truth rms
    model_pars = pack_model_pars(true_u, true_v,...
        true_A, true_B, true_lambda2_y);

    true_rms = get_validation_error(model_pars, pars, y_held,...
        d_held);
end
% plot(rms, '-*b')
% hold on
% plot(true_rms, '-+r')
bar([rms, 0], 'b')
hold on
bar([0, true_rms], 'r')
plot([noise, noise],'-g')
title(name);
legend('Learned',' True')
ylabel('RMS');
if savePlots
    saveas(gca(), ['rms', name, num2str(jj), '.pdf']);
%     hgexport(h, ['rms', name, num2str(jj), '.pdf']);
end

exitcode = 0;
end