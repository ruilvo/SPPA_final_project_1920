%% SPAA final assignment
% Rui Oliveira, 2020

% Exercise 2, verification of the first part

close all;
clear all;
clc;
format long

%% Parameteres
ns = 100000000; % Number of samples

% f = [1, 3, 3]; % v -> x impulse response (first case)
f=[3,3,1];  % v -> x impulse response (second case)
g = [1, -2, 1]; % v -> d impulse response

ks = -3:1:3; % Delays to try

%% Generate signals
v = randn(1, ns);
x = filter(f, 1, v);
d = filter(g, 1, v);

%% Simulation

%  x and x
% Store results
resultsxx = zeros(length(ks), 1);

for i = 1:length(ks)
    k = ks(i);

    if k >= 0
        % if k>=0 then delay d instead
        x1 = filter([zeros(1, abs(k)) 1], 1, x);
        % And do the math
        resultsxx(i) = dot(x, x1) / ns;
    else
        % delay x instead
        x1 = filter([zeros(1, abs(k)) 1], 1, x);
        resultsxx(i) = dot(x1, x) / ns;
    end % ks

end % ks

%  x and d
% Store results
resultsxd = zeros(length(ks), 1);

for i = 1:length(ks)
    k = ks(i);

    if k >= 0
        % if k>= then delay d instead
        d1 = filter([zeros(1, abs(k)) 1], 1, d);
        % And do the math
        resultsxd(i) = dot(x, d1) / ns;
    else
        % delay x instead
        x1 = filter([zeros(1, abs(k)) 1], 1, x);
        resultsxd(i) = dot(x1, d) / ns;
    end % ks

end % ks

resultsxx
resultsxd