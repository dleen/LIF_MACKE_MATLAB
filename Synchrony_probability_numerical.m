clear all;
clf

% Set the number of neurons used:
N = 100;

% Import the data to be used.
% Data file has the form:
% p(x_0) p(x_1) ... p(x_(N-1)) mu rho
%M = importdata('figure_2a_0.64.dat',' ');

M{1} = importdata('figure_2a_0.325.dat',' ');
M{2} = importdata('figure_2a_0.49.dat',' ');
M{3} = importdata('figure_2a_0.73.dat',' ');

CM = hsv(3);

% PP contains the probability distributions.
% Each row represents a single run at mean mu and corr rho.

for i=1:3
    PP = M{i}(:,1:N+1);
    mu = M{i}(:,N+2);
    rho = M{i}(:,N+3);

    % Plot data.
    semilogy(mean(PP),'color',CM(i,:));
    hold on
end

title('Numerical Simulation','fontsize',18)
fig2a_leg = legend('\rho = 0.05','\rho = 0.1','\rho = 0.25');
set(fig2a_leg,'FontSize',16);
xlabel('Population spike count i.e synchrony','fontsize',16)
ylabel('Probability','fontsize',16)
axis([0 N 1e-4 1]);