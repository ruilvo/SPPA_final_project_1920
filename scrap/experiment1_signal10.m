% SPAA final assignment
% Rui Oliveira, 2020

%% Configs
close all;
clear all;

% Simulation parameters
n = 10; % Order of NLMS
mu = 0.7; % Adaptation step
lambda = 0.99; % Error filter parameter
T = 10000; % Number of points
sigma = 0.01; % Error power
epsi = 1e-300; % Divison by zero protection

%% Initialization

% Create first signal to identify
x = randn(1, T);
d = filter([1, 2], [1, 0.5], x) + sigma * randn(1, T);

% For the sake of sanity, let's see the impulse response of the channel
h = filter([1, 2], [1, 0.5], [1 zeros(1, n)]);
% These are convenient places to make plots
T1 = round(0.1 * T);
T2 = round(0.5 * T);
T3 = round(0.6 * T);

% This is for storing the results
W = zeros(T + 1, n + 1); % Weights
xx = zeros(1, n + 1); % Tapped delay line signals
e = zeros(T, 1); % Error

%% Simulation

for t = 1:T

    % At the simulation mid-point, reduce step by 10
    if t == fix(T/2)
        mu = 0.1 * mu;
    end

    xx = [x(t) xx(1:n)]; % Tapped delay line signals
    y = W(t, :) * xx'; % Estimate the signal
    e(t) = d(t) - y; % Compute the error
    % Step into the next weights
    W(t + 1, :) = W(t, :) + 2 * mu * xx * e(t) / (epsi + xx * xx');
end % t = 1:T

% Filter the errors
xi = filter(1 - lambda, [1 -lambda], e.^2);

%% Plots

figure(1);

% Weights
subplot(2, 1, 1);
hold on;

% Plot the weights all in black
plot(0:T, W, 'k');

% Plot the theoretical weights
plot([0 T], [h; h]);

% Mark the positions for the weights
for i = 1:n
    text(1.01 * T, h(i), sprintf('w_%d', i - 1));
end

grid on
hold off

% Error
subplot(2, 1, 2);
hold on;

% Plot the filtered square error (dB)
plot(0:T - 1, 10 * log10(xi));

% Plot the average by sections
plot([T1-1 T2-1], 10 * log10(mean(xi(T1:T2))) * [1 1], 'k');
plot([T3-1 T-1], 10 * log10(mean(xi(T3:T))) * [1 1], 'k');
hold off
grid on
