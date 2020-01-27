%% partial solution of problem 5
%
clear all
close all

%% configuration
%
n=10;                     % transversal filter order
T=30000;                  % number of samples
f=[1,-1,2,-2,7,2,2,2,-1]; % channel impulse response
sigma=0.01;               % standard deviation of the channel noise
mu=0.1;                   % adaptation step size
M = 8; % PSK order

fname = "../results/ex5_1.mat";

%% Generate signals
data = randi([0 M-1],1,T);
u = pskmod(data,M,pi/M);
x=filter([1,-1,2,-2,7,2,2,2,-1],1,u)+ ...
       sigma*(sqrt(2)*randn(1,T)+1i*sqrt(2)*randn(1,T)); % received data
%%
W=zeros(T+1,n+1);
W(1,round(n/2))=1;
xx=zeros(1,n+1);
d=0*x;
e=zeros(T,1);
for t=1:T
  xx=[x(t) xx(1:n)]; % tapped delay line signals
  y=W(t,:)*xx';
  d(t) = y;
  e(t)= pskmod(pskdemod(y,M,pi/M),M,pi/M)-y;
  W(t+1,:)=W(t,:)+2*mu*xx*e(t)/(1e-300+xx*xx');
end

%% squared error
%
lambda=0.99; % forgetting factor
xi=filter(1-lambda,[1 -lambda],abs(e).^2);

%% cross correlation between u and d
%
r=zeros(1,2*n-1);
rn=100;
i=T-rn+1-n:T-n;
for k=-n:n
  r(n+1+k)=u(i+k)*d(i)';
end
[m,i]=max(abs(r));
delay=i-n-1

%% display
%

figure(1);
subplot(2,1,1);
plot(0:T,abs(W),'Linewidth',2);
title('Weights');
grid on
subplot(2,1,2);
plot(0:T-1,10*log10(xi),'Linewidth',2);
title('Squared error');
grid on

figure(2)
subplot(2,1,1);
plot(-n:n,abs(r),'.-','Linewidth',2);
title('Cross-correlation (modulus)');
grid on
subplot(2,1,2);
plot(-n:n,angle(r)/(2*pi),'.-','Linewidth',2);
title('Cross-correlation (normalized phase)');
grid on

figure(3); hold on;
scatter(real(d(1:100)),imag(d(1:100)));
scatter(real(d(end-100:end)),imag(d(end-100:end)));

save(fname,'W','xi',"d");
