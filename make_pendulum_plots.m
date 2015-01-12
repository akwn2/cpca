function failed = make_pendulum_plots(cc, truA, truB)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Default true matrices for pendula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
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
end

% Plot matrices A and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:3
    figure
    colormap winter
    caxis([-5,5])
%         suptitle(strcat('Case: ', num2str(ii), 'p_', cc));
    text(.75,1.25,strcat('Case: ', num2str(ii), 'p_', cc));
    for jj = 1:4
        load(strcat(num2str(ii),'p_',cc,'_out',num2str(jj),'.mat'));

        subplot(4, 4, 4 * jj - 3)
        imagesc(A)
        title('Learned A');
        xlabel('D');
        ylabel('M');

        subplot(4, 4, 4 * jj - 2)
        imagesc(truA(1:2 * ii, 1:jj))
        title('True A');
        xlabel('D');
        ylabel('M');

        subplot(4, 4, 4 * jj - 1)
        imagesc(B)
        title('Learned B');
        xlabel('D');
        ylabel('M');

        subplot(4, 4, 4 * jj)
        imagesc(truB(1:2 * ii, 1:jj))
        title('True B');
        xlabel('D');
        ylabel('M');
    end
    axis tight
end

% Plot error score (denoising)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
score = zeros(3,4);
for ii = 1:3
    for jj = 1:4
        load(strcat(num2str(ii),'p_c_out',num2str(jj),'.mat'));
        score(ii, jj) = run_stats{3}(end);
    end
    subplot(3,1,ii)
    bar(score(ii,:))
    title(strcat('Case: ', num2str(ii), 'p_', cc));
    xlabel('Latent space dimension');
    ylabel('Denoising error score');
end
failed = 0;
end

