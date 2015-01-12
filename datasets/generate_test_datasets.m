function retcode = generate_test_datasets(max_size)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 1
        max_size = 4;
    end

    for dd = 1:max_size
        [u_true, v_true, A_true, B_true, lambda2_y_true] = ...
            init_true_val(2 * dd, dd, 1, 0);

        y = generate_data_from_model(u_true, v_true, A_true, B_true, ...
            lambda2_y_true, 1300);

        save(strcat('./datasets/',num2str(dd),'p_m.mat'), ...
            'y','u_true', 'v_true','A_true','B_true','lambda2_y_true');
    end
    retcode = 0;
end

