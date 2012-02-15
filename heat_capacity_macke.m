clear all;
clf

% Set the number of neurons used:
N = 100;

% Import the data to be used.
% Data file has the form:
% p(x_0) p(x_1) ... p(x_(N-1)) mu rho
M = importdata('figure_2a_EIF_0.39.dat',' ');

% PP contains the probability distributions.
% Each row represents a single run at mean mu and corr rho.
PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

% Temperature
T = linspace(0.25,4,200);

% For a colorful plot.
CM = hsv(length(N));

for j=1:length(N)
    % Calculate the correct binomial coefficients.
    binom = zeros(N(j)+1,1);
    for i=0:N(j)
        binom(i+1) = nchoosek(N(j),i);
    end

    % Calculate the DG probability distribution.
    P = mean(PP)';

    % Calculate the heat capacity as a function of temperature.
    for i=1:length(T)
        % Function to raise the probability dist to the power 1/T.
        % Non-trivial.
        P_beta = prob_dist_power(P,1/T(i),binom);
                
        % Find non-zeros (because we will take the log of P_beta).
        ind = find(P_beta~=0);
        
        % Calculate the heat capacity. HC = var(log(p(x)))/n. 
        c(i) = sum(P_beta(ind).*(log2(P_beta(ind))-log2(binom(ind))).^2) -...
        sum(P_beta(ind).*(log2(P_beta(ind))-log2(binom(ind))))^2;
        c(i) = c(i)/N(j);
    end
    % Macke's plot has a log2 x axis.
    plot(log2(T),c,'color',CM(j,:));
    hold on
end

% Plot the line T=1.
plot(0,linspace(0,max(c),200),'r')
