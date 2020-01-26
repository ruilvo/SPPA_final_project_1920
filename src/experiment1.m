% SPAA final assignment
% Rui Oliveira, 2020

%% Configs
close all;
clear all;

%% Parameteres
T = 10000; % Number of points
lambda = 0.99;
nstry = [3, 5, 10];
mustry = [0.3, 0.6, 0.9];
sig = 0.01;

%% Simulations
% This is for storing the results
ds = zeros(3, 1, T);
W = zeros(3, length(nstry), length(mustry), T + 1, max(nstry)); % Weights
e = zeros(3, length(nstry), length(mustry), T, 1); % Error
ys = zeros(3, length(nstry), length(mustry), 1, T); % Computed signal

%% First signal
% Create signal
x = randn(1, T);
d = filter([1, 2], [1, 0.5], x) + sig * randn(1, T);
ds(1, :, :) = d;

for i = 1:length(nstry)
    n = nstry(i);

    for j = 1:length(mustry)
        mu = mustry(j);
        [ys(1, i, j, :, :), W(1, i, j, :, 1:n + 1), e(1, i, j, :, :)] = do_nlms(x, d, n, mu, lambda, T);

    end % length(mustry)

end % length(nstry)

%% Second signals
phase = 0.01 + 0.01 * floor(4 * (0:T - 1) / T);
x0 = randn(1, T);
x1 = sin(2 * pi * phase .* (0:T - 1));
x2 = cos(8 * pi * (0:T - 1) / T) >= 0;
x = x0 .* x2 + x1 .* (1 - x2);
d = filter([1, 2], [1, 0.5], x) + sig * randn(1, T);
ds(2, :, :) = d;

for i = 1:length(nstry)
    n = nstry(i);

    for j = 1:length(mustry)
        mu = mustry(j);
        [ys(2, i, j, :, :), W(2, i, j, :, 1:n + 1), e(2, i, j, :, :)] = do_nlms(x, d, n, mu, lambda, T);

    end % length(mustry)

end % length(nstry)

%% Third signal
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

d = d + sig * randn(1, T);
ds(3, :, :) = d;

for i = 1:length(nstry)
    n = nstry(i);

    for j = 1:length(mustry)
        mu = mustry(j);
        [ys(3, i, j, :, :), W(3, i, j, :, 1:n + 1), e(3, i, j, :, :)] = do_nlms(x, d, n, mu, lambda, T);

    end % length(mustry)

end % length(nstry)

%% Save data
save("../results/results1.mat", "nstry", "mustry", ...
    "W", "ys", "e", "ds");
