%test_do_e_step
% Runs 3 unit tests on the negative_log_joint function:
%   1.) Compare against a non-vectorised implementation
%   2.) Run with "true" mean field parameters

clear; clc;
fprintf('Testing do_e_step\n');
fprintf('--------------------------\n\n');

%% Set all parameters
%------------------------------
set_test_parameters;

%% Start test 1
%----------------------
eStepPasses = 10;
checkFq = false;
fq = -1E10;
fq_old = -2E10;
var_list = {'u','v','A','B','c','lambda2_y'};
[vararray, rows, cols] = packit(u, v, A, B, c, log(lambda2_y));

fprintf('\tTest 1: Compare with other implementation.\n');
[mc_loop, ms_loop, mc2_loop, msc_loop, ms2_loop, ...
 k1_loop, k2_loop, m1_loop, m2_loop, ...
 fq_loop, fq_old_loop] = do_e_step_loop(y, u, v, A, B, c, lambda2_y, ...
                                        mc, ms, mc2, ms2, msc, ...
                                        k1, k2, m1, m2, fq, fq_old,...
                                        eStepPasses, var_list, checkFq);

[mc_semi, ms_semi, mc2_semi, msc_semi, ms2_semi, ...
 k1_semi, k2_semi, m1_semi, m2_semi, ...
 fq_semi, fq_old_semi] = do_e_step_semi_vec(y, u, v, A, B, c, lambda2_y, ...
                                 mc, ms, mc2, ms2, msc, ...
                                 k1, k2, m1, m2, fq, fq_old,...
                                 eStepPasses, var_list, checkFq);

[mc_vec, ms_vec, mc2_vec, msc_vec, ms2_vec, ...
 k1_vec, k2_vec, m1_vec, m2_vec, ...
 fq_vec, fq_old_vec] = do_e_step(y, u, v, A, B, c, lambda2_y, ...
                                 mc, ms, mc2, ms2, msc, ...
                                 k1, k2, m1, m2, fq, fq_old,...
                                 eStepPasses, var_list, checkFq);

%% Plots of test 1
% reshaping everything to be able to plot nicely
mc_loop = reshape(mc_loop, 1, D * N);
ms_loop = reshape(ms_loop, 1, D * N);
mc2_loop = reshape(mc2_loop, 1, D * N);
ms2_loop = reshape(ms2_loop, 1, D * N);
msc_loop = reshape(msc_loop, 1, D * N);
k1_loop = reshape(k1_loop, 1, D * N);
k2_loop = reshape(k2_loop, 1, D * N);
m1_loop = reshape(m1_loop, 1, D * N);
m2_loop = reshape(m2_loop, 1, D * N);

mc_semi = reshape(mc_semi, 1, D * N);
ms_semi = reshape(ms_semi, 1, D * N);
mc2_semi = reshape(mc2_semi, 1, D * N);
ms2_semi = reshape(ms2_semi, 1, D * N);
msc_semi = reshape(msc_semi, 1, D * N);
k1_semi = reshape(k1_semi, 1, D * N);
k2_semi = reshape(k2_semi, 1, D * N);
m1_semi = reshape(m1_semi, 1, D * N);
m2_semi = reshape(m2_semi, 1, D * N);


mc_vec = reshape(mc_vec, 1, D * N);
ms_vec = reshape(ms_vec, 1, D * N);
mc2_vec = reshape(mc2_vec, 1, D * N);
ms2_vec = reshape(ms2_vec, 1, D * N);
msc_vec = reshape(msc_vec, 1, D * N);
k1_vec = reshape(k1_vec, 1, D * N);
k2_vec = reshape(k2_vec, 1, D * N);
m1_vec = reshape(m1_vec, 1, D * N);
m2_vec = reshape(m2_vec, 1, D * N);

% First figure: expectations of trigonometric functions under the mean
% field

figure
subplot(2,5,1)
plot(mc_vec, mc_loop,'o');
hold all
plot(mc_vec, mc_semi,'o');
plot(-1:1,-1:1);
subplot(2,5,2)
plot(ms_vec, ms_loop,'o');
hold all
plot(ms_vec, ms_semi,'o');
plot(-1:1,-1:1);
subplot(2,5,3)
plot(mc2_vec, mc2_loop,'o');
hold all
plot(mc2_vec, mc2_semi,'o');
plot(-1:1,-1:1);
subplot(2,5,4)
plot(ms2_vec, ms2_loop,'o');
hold all
plot(ms2_vec, ms2_semi,'o');
plot(-1:1,-1:1);
subplot(2,5,5)
plot(msc_vec, msc_loop,'o');
hold all
plot(msc_vec, msc_semi,'o');
plot(-1:1,-1:1);

% Differences
subplot(2,5,6)
plot((mc_vec - mc_loop) ./ abs(mc_vec));
hold all
plot((mc_vec - mc_semi) ./ abs(mc_vec));
subplot(2,5,7)
plot((ms_vec - ms_loop) ./ abs(ms_vec));
hold all
plot((ms_vec - ms_semi) ./ abs(ms_vec));
subplot(2,5,8)
plot((mc2_vec - mc2_loop) ./ abs(mc2_vec));
hold all
plot((mc2_vec - mc2_semi) ./ abs(mc2_vec));
subplot(2,5,9)
plot((ms2_vec - ms2_loop) ./ abs(ms2_vec));
hold all
plot((ms2_vec - ms2_semi) ./ abs(ms2_vec));
subplot(2,5,10)
plot((msc_vec - msc_loop) ./ abs(msc_vec));
hold all
plot((msc_vec - msc_semi) ./ abs(msc_vec));

% Second plot: Mean field parameters
%------------------------

figure
subplot(2,4,1)
plot(k1_vec, k1_loop,'o');
hold all
plot(k1_vec, k1_semi,'o');
plot(0:max(max([k1_loop, k1_vec, k1_semi])),...
    0:max(max([k1_loop, k1_vec, k1_semi])));
subplot(2,4,2)
plot(k2_vec, k2_loop,'o');
hold all
plot(k2_vec, k2_semi,'o');
plot(0:max(max([k2_loop, k2_vec, k2_semi])),...
    0:max(max([k2_loop, k2_vec, k2_semi])));
subplot(2,4,3)
plot(m1_vec, m1_loop,'o');
hold all
plot(m1_vec, m1_vec,'o');
plot(-pi:pi,-pi:pi);
subplot(2,4,4)
plot(m2_vec, m2_loop,'o');
hold all
plot(m2_vec, m2_semi,'o');
plot(-pi:pi,-pi:pi);

subplot(2,4,5)
plot((k1_vec - k1_loop)./abs(k1_vec));
hold all
plot((k1_vec - k1_semi)./abs(k1_vec));
subplot(2,4,6)
plot((k2_vec - k2_loop)./abs(k2_vec));
hold all
plot((k2_vec - k2_semi)./abs(k2_vec));
subplot(2,4,7)
plot((m1_vec - m1_loop)./abs(m1_vec));
hold all
plot((m1_vec - m1_semi)./abs(m1_vec));
subplot(2,4,8)
plot((m2_vec - m2_loop)./abs(m2_vec));
hold all
plot((m2_vec - m2_semi)./abs(m2_vec));