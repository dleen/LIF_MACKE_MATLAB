clear all;

% Set the number of neurons used:
N = 100;

% Import the data to be used.
% Data file has the form:
% p(x_0) p(x_1) ... p(x_(N-1)) mu rho
M = importdata('figure_2a_0.73.dat',' ');

% PP contains the probability distributions.
% Each row represents a single run at mean mu and corr rho.
PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

% Plot data.
semilogy(mean(PP),'b');
hold on
axis([0 N 1e-4 1]);