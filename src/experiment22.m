%% partial solution of the first part of problem 2
%%
close all;
clear all;
%% configuration
%%
N = 60; % largest possible filter order
% f = [1, 3, 3]; % v -> x impulse response (first case)
f = [3, 3, 1]; % v -> x impulse response (second case)
g = [1, -2, 1]; % v -> d impulse response
fname = "../results/results22_2.mat"; % adjust accordingly

%% theoretical inner products
%%
rxx = conv(f, f(end:-1:1)); % complete auto-correlation
rxd = conv(g, f(end:-1:1)); % complete cross-correlation
rdd = sum(g.^2); % variance of d

% get the parts that will be needed later
rxx = [rxx(length(f):end), zeros(1, N + 1 - length(f))];
rxd = [rxd(length(f):end), zeros(1, N + 1 - length(f))];

%% optimal solution for a specific order of the transversal filter
%%
close all
EDb = zeros(1, N + 1);
W = zeros(N + 1, N + 1);

for n = 1:N + 1
    R = toeplitz(rxx(1:n)); % R matrix
    p = rxd(1:n)'; % p vector
    w = pinv(R) * p; % optimal weights
    E = rdd - w' * p; % mimimum squared error
    EDb(n) = 10 * log10(E);
    W(n, 1:n) = w';
end

figure(1); hold on;
plot(1:N + 1, EDb, '.-');
grid on;
figure(2); hold on;

for n = 1:N + 1
    plot([n:N+1].', W(n:end, n), '.-');
end

grid on;

%% quality control (optional)
%%
if 0
    ns = 1000000; % number of samples
    v = randn(1, ns)';
    x = filter(f, 1, v);
    d = filter(g, 1, v);
    X = toeplitz(x, [x(1) zeros(1, n)]);
    X' * X / ns - R
    X' * d / ns - p
    pinv(X) * d - w
end

%% estimation of the optimal IIR filter
%%
if 0
    W = toeplitz(w, [w(1) 0 0]);
    svd(W(10:end, :))
    A = W(4:end, 1:2);
    b = W(4:end, 3);
    d = pinv(A) * b;
    d = [d; -1];
    W * d
end

%% Save data
save(fname, "W", "EDb", "N");
