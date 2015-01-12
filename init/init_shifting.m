function [u, v, A, B, lambda2_y] = init_shifting(y, u, v, A, B, lambda2_y)
%INIT_SUPERVISED Initializes the problem using several supervised problems
%   This initialization procedure is based on displacing each pendulum's
%   base to the origin and learning a representation for it disconsidering
%   the cross effects from other pendula. To do this we assume that y is
%   ordered in such a way that y(1:2,:) is the first pendulum, y(3:4,:) is
%   the second pendulum and so on and that y(2 * ii - 1:2 * ii, :) is
%   dependent only on the previous pendulum position, that is 
%   y(2 * ii - 3 : 2 * ii - 2, :) for ii = 2, 3, ..., n_dim. The priors are
%   then calculated by fitting the noisy angles to a von Mises.
%   The idea can be summarised loosely as a Grahm-Schmidt procedure on the
%   pendula.

% Fix the case where the dimensionality of the hidden space is larger than
% the observed space (those will be set to zero)
[m_dim, d_pts] = size(y);
n_dim = size(A, 2);

if m_dim < 2 * n_dim
    n_dim = floor(m_dim / 2);
end

% Pre allocation - prior scaling must not be zero or causes problems when
% computing values
p_scale = 1.E-6;
u = p_scale .* ones(size(u));
v = p_scale .* ones(size(v));
A = zeros(size(A));
B = zeros(size(B));
phi = zeros(size(A,2), d_pts);

subA = zeros(m_dim, n_dim);
subB = zeros(m_dim, n_dim);
theta = zeros(n_dim, d_pts);
sigma2 = zeros(2 * n_dim, 1);

for ii = 1:n_dim
    % Displacement to origin of complex plane
    if ii == 1
        z = y(1,:) + 1.i * y(2,:);
    else
        z = (y(2 * ii - 1, :) - y(2 * ii - 3,:)) + ...
                1.i * (y(2 * ii,:) - y(2 * ii - 2,:));
    end
    % Noisy and angle and absolute value
    theta(ii, :) = angle(z);
    subA(2 * ii - 1:2:m_dim, ii) = mean(abs(z));
    subB(2 * ii:2:m_dim, ii) = mean(abs(z));
    sigma2(ii) = var(abs(z));
end

% Fit von Mises to noisy angles using estimators in the form of maximum
% likelihood, in particular the estimator for kappa from Dobson (1978)
S = sin(theta);
C = cos(theta);
R = (mean(C, 2).^2 + mean(S, 2).^2).^0.5;
mu = acos(mean(C, 2));
kappa = (1.28 - 0.53 * R.^2) .* tan(0.5 .* pi .* R);

% Output all model parameters
u(1:n_dim) = kappa .* cos(mu);
v(1:n_dim) = kappa .* sin(mu);

A(:, 1:n_dim) = subA;
B(:, 1:n_dim) = subB;

% if max(size(lambda2_y)) == 1
%     phi(1:n_dim, :) = theta;
%     aux = sqrt(sum(sum((y - A * cos(phi) + B * sin(phi)).^2)));
%     lambda2_y = (m_dim * d_pts - 1) / aux;
% else
    lambda2_y = 1 ./ max(sigma2);
% end

end