%test_negative_log_joint
% Runs 3 unit tests on the negative_log_joint function:
%   1.) Compare against a non-vectorised implementation
%   2.) Run checkgrad
%   3.) Run an optimisation

clear; clc;
fprintf('Testing negative_log_joint\n');
fprintf('--------------------------\n\n');

%% Set all parameters
%------------------------------
set_test_parameters;

%% Start test 1
%----------------------
fprintf('\tTest 1: Compare with other implementation.\n');
loop_log_p = negative_log_joint_loop(y, mc, ms, mc2, ms2,...
                                     msc, u, v, A, B, c, lambda2_y);
fprintf('\t\tLoop implementation of -<log p> = %1.4e\n', loop_log_p);

var_list = {'u','v','A','B','c','lambda2_y'};
[vararray, rows, cols] = packit(u, v, A, B, c, log(lambda2_y));
        
semi_log_p = negative_log_joint_semi_vec(vararray, var_list, rows, cols, ...
    y, mc, ms, mc2, ms2, msc, ones(size(u)), ones(size(v)),...
    zeros(size(A)), zeros(size(B)), zeros(size(c)), ones(size(lambda2_y)));
fprintf('\t\tSemi-vectorised implementation of -<log p> = %1.4e\n', semi_log_p);

full_log_p = negative_log_joint(vararray, var_list, rows, cols, ...
    y, mc, ms, mc2, ms2, msc, ones(size(u)), ones(size(v)),...
    zeros(size(A)), zeros(size(B)), zeros(size(c)), ones(size(lambda2_y)));

fprintf('\t\tFully vectorised implementation of -<log p> = %1.4e\n', full_log_p);

passed1 = norm(loop_log_p - full_log_p) < tol && norm(semi_log_p - full_log_p) < tol;
if ~passed1
    passed1 = false;
    fprintf('\t\t-------------------------\n');
    fprintf('\t\t!!! TEST 1 NOT PASSED !!!\n');
    fprintf('\t\t-------------------------\n');
    fprintf('\t\t||loop - vectorised|| = %1.4e\n', norm(loop_log_p - full_log_p))
    fprintf('\t\t||semi - vectorised|| = %1.4e\n', norm(semi_log_p - full_log_p))
    fprintf('\t\t||semi - loop|| = %1.4e\n', norm(semi_log_p - loop_log_p))
else
    fprintf('\t\tTest 1: OK.\n');
end

%% Start test 2
%----------------------
fprintf('\tTest 2: Run checkgrad on gradients.\n');
e = 1.E-8;
dtol = 1.E-5;
gradients = 'negative_log_joint';
fprintf(['\t\t\tLog-joint function used: ', gradients, '\n']);

var_list = {'u','v','A','B','c', 'lambda2_y'};
passed2 = true;
for ii = 1:length(var_list)
    fprintf(['\t\t\t Gradients of ', var_list{ii}, '\n']);
    
    [vararray, rows, cols] = packit(eval(var_list{ii}));
    
    d = checkgrad(gradients, vararray', e,...
        var_list(ii), rows, cols, ...
        y, mc, ms, mc2, ms2, msc, ones(size(u)), ones(size(v)),...
        zeros(size(A)), zeros(size(B)), zeros(size(c)), ones(size(lambda2_y)));
    
    fprintf('\t\tRelative difference: %1.4e\n', d);
    passed2 = passed2 && (d < dtol);
end

if ~passed2
    fprintf('\t\t-------------------------\n');
    fprintf('\t\t!!! TEST 2 NOT PASSED !!!\n');
    fprintf('\t\t-------------------------\n');
else
    fprintf('\t\tTest 2: OK.\n');
end


%% Test 3
% Optimisation step
fprintf('\tTest 3: Optimisation.\n');

gradients = 'negative_log_joint';
fprintf(['\t\t\tLog-joint function used: ', gradients, '\n']);

var_list = {'u','v','A','B','c', 'lambda2_y'};
% passed3 = true;
mStepMaxIter = 1E3;
for ii = 1:length(var_list)
    fprintf(['\t\t\t Optimising over ', var_list{ii}, '\n']);
    
    [vararray_true, rows, cols] = packit(eval(var_list{ii}));
    vararray_noisy = vararray_true + ...
                    0.5 * max(vararray_true) * rand(size(vararray_true));
    [vararray, fX, iter] = minimize(vararray_noisy,...
                                    gradients, - mStepMaxIter, ...
                                    var_list(ii), rows, cols,...
                                    y, mc, ms, mc2, ms2, msc, ...
                                    u, v, A, B, c, lambda2_y);
    
    fprintf('\t\t\tNumber of iterations: %d.\n', length(fX));
    fprintf('\t\t\tNumber of line searches: %d.\n', iter);
    d = norm(vararray - vararray_true);
    fprintf('\tStart\tFinal\tTrue\n');
    disp([vararray_noisy', vararray', vararray_true']);
    fprintf('\t\t\tNorm difference: %1.4e\n', d);

%     passed3 = passed3 && (d < dtol);
end

% if ~passed3
%     fprintf('\t\t-------------------------\n');
%     fprintf('\t\t!!! TEST 3 NOT PASSED !!!\n');
%     fprintf('\t\t-------------------------\n');
% else
%     fprintf('\t\tTest 4: OK.\n');
% end
                        
                        
%% Test summary
fprintf('-------------------------\n');
fprintf('      TESTS SUMMARY      \n');
fprintf('-------------------------\n');

if passed1
    fprintf('Test 1: OK.\n');
else
    fprintf('Test 1: Failed.\n');
end

if passed2
    fprintf('Test 2: OK.\n');
else
    fprintf('Test 2: Failed.\n');
end