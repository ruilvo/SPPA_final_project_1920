close all;
clear all;

n = 2; % filter order
miu = 0.01; % step size 

T=10000;
sigma=0.01;
x=randn(1,T);
d=filter([1,2],[1,0.5],x)+sigma*randn(1,T);

% Transversaal filter order n
%          ------        ------
%  x ---- | z^-1 | ---- | z^-1 | --- ... ---
%     |    ------    |   ------           |
%     |              |                    |
%     x0             x1                  xn

% Theoretical impulse response
h = filter([1,2],[1,0.5], [1 zeros(1,n)]);

% time vertically 
% weights horizontally
W = zeros(T+1, n+1);

xx = zeros(1, n+1);
e = zeros(T,1);

for t=1:T
    if t==floor(T/2)
        miu = 0.1*miu;
    end
    xx =[x(t) xx(1:n)]; % tapped delay line signals
    y = W(t,:)*xx';
    e(t) = d(t)-y;
    W(t+1,:) = W(t,:)+2*miu*xx*e(t) / (1e-300 + xx*xx');
end

lambda = 0.99; % Forget factor
% winfilter = filter(1-lambda

figure(1); 
subplot(2,1,1);
plot(0:T,W);
hold on; grid on;
for i=1:n+1
    yline(h(i));
end
hold off;
subplot(2,1,2);

