clear all;

N = 100;

tic
M = importdata('figure_2a_0.2.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

semilogy(mean(PP),'r');