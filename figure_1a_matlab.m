clear all;
clf

% Number of neurons.
N = 5;

% Number of samples to average over.
avnum = 200;

% Import data for use in plots.
M{1} = importdata('figure_1a_0.1.dat',' ');
M{2} = importdata('figure_1a_0.3.dat',' ');
M{3} = importdata('figure_1a_0.5.dat',' ');

% For colorful plots.
CM = hsv(3);

% Loop for each of the lines in the plot.
for j=1:3
    % The means.
    mu  = M{j}(:,N+2);
    % The correlation coefficients.
    rho = M{j}(:,N+3);

    mmu  = zeros(length(M{j})/avnum,1);
    rrho = zeros(length(M{j})/avnum,1);

    % Find any NaNs.
    rho(isnan(rho)) = 0;

    % Sort the means.
    [mu_sorted,IX] = sort(mu);
    rho_sorted = rho(IX);

    % Average the means and correlations.
    for i=1:length(M{j})/avnum
        mmu(i) = mean(mu_sorted((avnum*(i-1)+1):(avnum*i)));
        rrho(i) = mean(rho_sorted((avnum*(i-1)+1):(avnum*i)));
    end

    plot(mmu,rrho,'color',CM(j,:));
    hold on
end

xlabel('Mean \mu','fontsize',16)
ylabel('Correlation \rho','fontsize',16)
fig1a_leg = legend('\lambda = 0.1','\lambda = 0.3','\lambda = 0.5');
set(fig1a_leg,'fontsize',16)
axis([0 0.5 0 0.2])