% A study of the Generalised von Mises.

% % Analyse GvM to conditioned Gaussian
% S2X = 1.0;
% S2Y = 1.0;
% RHO = 0.5;
% C = [S2X, RHO; RHO, S2Y];
% CInv = C.^-1;
% 
% r = 1.0;
% theta = -pi:0.01:pi;
% X = r * cos(theta);
% Y = r * sin(theta);
% 
% rm = 0.8;
% thetam = pi/2;
% MX = rm * cos(thetam);
% MY = rm * sin(thetam);
% 
% GvM = exp(- (X - MX) .^ 2 * CInv(1,1) + ...
%     (Y - MY) .^ 2 * CInv(2,2) + ...
%     (X - MX) .* (Y - MY) .* (C(1,2) + C(2,1))) / sqrt((2 * pi)^2 * det(C));
% 
% plot3(X, Y, GvM);
% grid on;

% Analyse our approximation
k1 = 60.0;
k2 = 30.0;
m1 = 0.;
m2 = 1.5;

theta = -pi:0.01:pi;

u_gvm = exp(k1 * cos(theta - m1) + k2 * cos(2*(theta - m2)));
z_gvm = sum(u_gvm);

gvm = u_gvm ./ z_gvm;

figure
plot(theta, gvm);

mode = gvm_modes(k1, k2, m1, m2);

hold on
stem(mode);