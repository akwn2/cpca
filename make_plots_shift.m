function exitcode = make_plots_shift(cc, save_plots, truA, truB, truVar)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    save_plots = 0;
end
% Default true matrices for pendula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    truA = [1, 0, 0, 0;
            0, 0, 0, 0;
            1, 1, 0, 0;
            0, 0, 0, 0;
            1, 1, 1, 0;
            0, 0, 0, 0;
            1, 1, 1, 1;
            0, 0, 0, 0];

    truB = [0, 0, 0, 0;
            1, 0, 0, 0;
            0, 0, 0, 0;
            1, 1, 0, 0;
            0, 0, 0, 0;
            1, 1, 1, 0;
            0, 0, 0, 0;
            1, 1, 1, 1];

    truA = truA .* 5;
    truB = truB .* 5;
    truVar = 1;
end

% Plot matrices A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:3
    figure
    %colormap summer
    clim = [-5, 15];
    
    caxis([-5,5])
%         suptitle(strcat('Case: Shift  ', num2str(ii), 'p_', cc));
    text(.75, 1.25, strcat('Shift Case: Shift  ', num2str(ii), 'p_', cc));
    ax = zeros(16,1);
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'shift.mat'));

        ax(4 * jj - 3) = subplot(4, 4, 4 * jj - 3);
        imagesc(A, clim)
        title('Learned A');
        xlabel('D');
        ylabel('M');

        ax(4 * jj - 2) = subplot(4, 4, 4 * jj - 2);
        imagesc(truA(1:2 * ii, 1:jj), clim)
        title('True A');
        xlabel('D');
        ylabel('M');

        ax(4 * jj - 1) = subplot(4, 4, 4 * jj - 1);
        imagesc(B, clim)
        title('Learned B');
        xlabel('D');
        ylabel('M');

        ax(4 * jj) = subplot(4, 4, 4 * jj);
        imagesc(truB(1:2 * ii, 1:jj), clim)
        title('True B');
        xlabel('D');
        ylabel('M');
    end
    axis tight
    gcf();
    h = colorbar;
    set(h, 'Position', [.90 .125 .0525 .8])
    
    for jj=1:16
          pos=get(ax(jj), 'Position');
          set(ax(jj), 'Position', [pos(1) pos(2) 0.75 * pos(3) pos(4)]);
    end
    
    if save_plots
        saveas(gca(), strcat('./figs/matrices_case', num2str(ii), 'p_', cc, 'shift.fig'));
    end
end



% Plot variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
score = ones(4,1);
for ii = 1:3
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'shift.mat'));
        score(jj) = 1 / lambda2_y;
    end
    subplot(3, 1, ii)
    bar(score)
    title(strcat('Case: Shift  ', num2str(ii), ' p ', cc));
    ylim([0, max([truVar, score']) + 1.]);
    ylabel('Variance');
    xlabel('Latent space size');
    hold on
    plot(0:5, truVar .* ones(6,1));
end
if save_plots
    saveas(gca(), strcat('./figs/denoising_p_', cc, 'shift.fig'));
end



% Plot error score (denoising)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
score = zeros(3,4);
for ii = 1:3
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'shift.mat'));
        
        [m_dim, n_dim] = size(A);
        d_held = size(y_held, 2);
    
        model_pars = pack_model_pars(u, v, A, B, lambda2_y);
    
        pars = pack_pars(y_held, rand(n_dim, d_held), rand(n_dim, d_held), ...
            rand(n_dim, d_held), rand(n_dim, d_held), rand(n_dim, d_held), ...
            n_dim, m_dim, d_held);
        
        score(ii, jj) = get_validation_error(model_pars, pars, y_held, d_held);
    end
    subplot(3,1,ii)
    bar(score(ii,:))
    title(strcat('Case: Shift  ', num2str(ii), 'p ', cc));
    xlabel('Latent space dimension');
    ylabel('RMS');
end
if save_plots
    saveas(gca(), strcat('./figs/denoising_p_', cc, 'shift.fig'));
end



% Plot model versus error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ii = 1:3
%     for jj = 1:4
%         load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'shift.mat'));
%         [m_dim, d_pts] = size(y);
%         
%         [~, y_model] = generate_data_from_model(u, v, A, B,...
%                                                 lambda2_y, d_pts);
%         figure
%         text(.75, 1.25, strcat('Case: Shift  ', num2str(ii), 'p_', cc,...
%             num2str(jj)));
%         
%         for kk = 1:m_dim
%             subplot(2,m_dim / 2, kk)
%             plot(y, y_model)
%             xlabel(strcat('True value for y(', num2str(kk),',:)'))
%             ylabel(strcat('Learned value for y(', num2str(kk),',:)'))
%         end
%         if save_plots
%             saveas(gca(), strcat('./figs/cross_case', num2str(ii), ...
%                 'p_', cc, num2str(jj), 'shift.fig'));
%         end
%     end
% end

exitcode = 0;
end