% SPAA ex 5
% Rui Oliveira
% map-tele 2020

clc;
clear all;
close all;

%% Configuration
T = 10000;
omega = 0.2*pi;
sigm = 0.01;

%% Signals
tvec = 0:T+1;
x = cos(omega*tvec)+sigm*rand(1,T+2);

d1 = x(2:end); % For the delay by one
% Or the second case
A = 2;
d2 = sin(omega*tvec).*(A*(tvec<T/2)-A*(tvec>=T/2)); 

% And remove the ones I put extra
x = x(1:end-1); 
tvec = tvec(1:end-1);
d2 = d2(1:end-1);

%% Simulation

% Select here which one
d = d1;

for mu = [0.1 0.8]
    for n = [1 2 3]
        fname = sprintf("../results/results_ex4_d1_n%d_mu%.2f.mat",n,mu);

        % This is for storing the results
        W = zeros(T + 2, n); % Weights
        xx = zeros(1, n); % Tapped delay line signals
        e = zeros(T+1, 1); % Error
        ys = zeros(1, T+1); % Computed signal

        epsi = 1e-300; % Divison by zero protection

        for t = 1:T+1
            xx = [x(t) xx(1:n-1)]; % Tapped delay line signals

            y = W(t, :) * xx'; % Estimate the signal
            ys(t) = y;
            e(t) = d(t) - y; % Compute the error

            % Step into the next weights
            W(t + 1, :) = W(t, :) + 2 * mu * xx * e(t) / (epsi + xx * xx');

        end % for t

        lambda = 0.95;
        err = filter(1 - lambda, [1 -lambda], abs(e).^2);

        save(fname, "W","err","n","mu","d","ys");
    end
end

% Select here which one
d = d2;

for mu = [0.1 0.8]
    for n = [1 2 3]
        fname = sprintf("../results/results_ex4_d2_n%d_mu%.2f.mat",n,mu);

        % This is for storing the results
        W = zeros(T + 2, n); % Weights
        xx = zeros(1, n); % Tapped delay line signals
        e = zeros(T+1, 1); % Error
        ys = zeros(1, T+1); % Computed signal

        epsi = 1e-300; % Divison by zero protection

        for t = 1:T+1
            xx = [x(t) xx(1:n-1)]; % Tapped delay line signals

            y = W(t, :) * xx'; % Estimate the signal
            ys(t) = y;
            e(t) = d(t) - y; % Compute the error

            % Step into the next weights
            W(t + 1, :) = W(t, :) + 2 * mu * xx * e(t) / (epsi + xx * xx');

        end % for t

        lambda = 0.95;
        err = filter(1 - lambda, [1 -lambda], abs(e).^2);

        save(fname, "W","err","n","mu","d","ys");
    end
end