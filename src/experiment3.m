%% SPAA Final assignment, exercise 3
% NLMS algorithm for a transversal filter
% inputs:
% x in the input signal
% d is the desired signal
% n is the transversal filter order
% mu is the convergence factor
% outputs:
% W are the weights
% y is the output
% e is the error

clc;
close all;
clear all;

%% Configuration
ns = 2000000; % number of samples
f = [1, 3, 3]; % v -> x impulse response
% f = [3, 3, 1]; % v -> x impulse response, 2nd signal
g = [1, -2, 1]; % v -> d impulse response
n = 40; %
mu = 0.1; %
epsi = 1e-100;
fname = "../results/results3_1.mat"; % adjust accordingly

%% Signals

tvec = 1:ns;
v = randn(ns, 1);

x = filter(f, 1, v);
d = filter(g, 1, v);

%% Simulation

mufac = ns / 1000;

T = length(x);
W = zeros(n + 1, T + 1);
xx = zeros(n + 1, 1);
y = zeros(size(x));
% d=zeros(size(x));
e = zeros(T, 1);
us = zeros(T, 1);

for t = 1:T
    mu = 0.1^(log10(10 + t / mufac));
    us(t) = mu;
    xx = [x(t); xx(1:n)];
    y(t) = W(:, t)' * xx;
    e(t) = d(t) - y(t);
    W(:, t + 1) = W(:, t) + 2 * mu / (epsi + xx' * xx) * e(t) * xx;
end

%% Plots
% Filter the error
lam = 0.99;
errav = filter([1 - lam, 0], [1, -lam], e.^2);

figure(1)
plot(10 * log10(errav));

figure(2); hold on;

for i = 1:n
    plot(W(i, :))
end

mu(end)
format long
W(:,end)

%% Save data
% save(fname, "W", "errav", "ns", "n");
