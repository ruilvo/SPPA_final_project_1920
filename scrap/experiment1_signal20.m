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
d = zeros(1, T);
U = zeros(4, T);
xx = [0, 0, 0, 0];

for t = 1:T
    xx = [x(t), xx(1:3)];
    U(:, t) = [0.2 + 0.3 * t / T; ...
                0.7 * cos(4 * pi * t / T); ...
                0.3 * sign(sin(10 * pi * t / T)) - 0.1; ...
                -0.2];
    d(t) = xx * U(:, t);
end

d = d + sigma * randn(1, T);

% These are convenient places to make plots
T1 = round(0.1 * T);
T2 = round(0.5 * T);
T3 = round(0.6 * T);

% This is for storing the results
W = zeros(T + 1, n + 1); % Weights
xx = zeros(1, n + 1); % Tapped delay line signals
e = zeros(T, 1); % Error
ys = zeros(1, T); % Computed signal

%% Simulation

for t = 1:T

    % At the simulation mid-point, reduce step by 10
    if t == T2
        mu = 0.1 * mu;
    end

    xx = [x(t) xx(1:n)]; % Tapped delay line signals
    y = W(t, :) * xx'; % Estimate the signal
    ys(t) = y;
    e(t) = d(t) - y; % Compute the error
    % Step into the next weights
    W(t + 1, :) = W(t, :) + 2 * mu * xx * e(t) / (epsi + xx * xx');
end % t = 1:T

% Filter the errors
xi = filter(1 - lambda, [1 -lambda], e.^2);

%% Plots

figure(1);

% Signals
subplot(3, 1, 1);
hold on;
plot(d);
plot(ys);
plot(d - ys);

legend("Signal",'Recovered',"Difference")

% Weights
subplot(3, 1, 2);
hold on;

% Plot the weights
plot(0:T, W);

% Mark the positions for the weights
for i = 1:n
    text(1.01 * T, W(end, i), sprintf('w_%d', i - 1));
end

grid on
hold off

% Error
subplot(3, 1, 3);
hold on;

% Plot the filtered square error (dB)
plot(0:T - 1, 10 * log10(xi));

hold off
grid on
