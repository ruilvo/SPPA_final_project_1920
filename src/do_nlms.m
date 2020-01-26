function [ys, W, err] = do_nlms(x, d, n, mu, lambda, T)
    % Execute the NLMS
    % d - signal to compute
    % n - order
    % miu - step
    % lambda - Filter for the error squared
    % T - Number of steps

    epsi = 1e-300; % Divison by zero protection

    % This is for storing the results
    W = zeros(T + 1, n + 1); % Weights
    xx = zeros(1, n + 1); % Tapped delay line signals
    e = zeros(T, 1); % Error
    ys = zeros(1, T); % Computed signal

    % x = randn(1, T);

    for t = 1:T

        % At the simulation mid-point, reduce step by 10
        if t == fix(T / 2)
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
    err = filter(1 - lambda, [1 -lambda], e.^2);
end
