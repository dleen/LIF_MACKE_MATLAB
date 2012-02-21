clear all;
tic

CM = hsv(4);

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 2;

% Number of neurons.
N = 5;
% Number of repeated trials at fixed mu, rho.
avnum = 20;

% Load data from a file. 
% The file.dat has the form:
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
% We keep mu and rho fixed over trials. We have multiple
% trials so that we can get an average and smooth out
% some of the inherent noise.
M{1} = importdata('figure_1b_LIF_0.2.dat',' ');
M{2} = importdata('figure_1b_-60.5.dat',' ');
M{3} = importdata('figure_1b_-61.2.dat',' ');
M{4} = importdata('figure_1b_-62.dat',' ');

% Generate a list of binomial coefficients.
binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

% List the possible states, using the symmetry argument.
% Symmetry argument can be found in Macke 2011.
% Essentially says that the parameters in the Ising model h_i, J_ij
% using symmetry conditions can be reduced to just 2 parameters. The 
% symmetry condition means that we talk about population spike 
% counts rather than particular probabilities. i.e P(101101) becomes 
% just P(4) as 4 spikes occurred. 
% Here we list the possible states which are just 1...N and at the
% same time we list their squares as needed in the Ising model.
states = [(0:N)',(0:N)'.^2];

% An initial guess:
param_list_init = [0;0];

for j=1:4
    % PP is a matrix where each row is a prob dist for a given
    % value of mu and rho.
    PP  = M{j}(:,1:N+1);
    mu  = M{j}(:,N+2);
    rho = M{j}(:,N+3);

    DKL = zeros(length(M{j}),1);

    for k=1:length(M{j})
        P = (PP(k,:))';
        % Calculate the means and moments i.e. the quantities
        % sum(k*P_k) and sum(k^2*P_k).
        mean_feature = P'*states;

        if alg_choice == 1
            % MATLAB Optimization Toolbox minimization function.
            % Uses Optimization toolbox unconstrained minimization function.
            options = optimset('GradObj','on','LargeScale','on',...
            'Display','final-detailed',...
            'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
            % We set GradObj = on as we supply the gradient.
            % We set LargeScale = on as that is the algorithm that uses the
            % user supplied gradient.
            % Display just shows some accuracy output.
            % The final options are tolerances etc.

            % The minimization function
            [param_list,fval,output] =...
                fminunc(@(x)neg_log_like_binom...
                (x,states,mean_feature,P,binom),param_list_init,options);
        elseif alg_choice == 2
            % Eric's minimization function, based on conjugate gradient.
            param_list = minimize(param_list_init,...
                'neg_log_like_binom',200,...
                states,mean_feature,P,binom);
        end

        % Use the parameters found via minimization to calculate the
        % probability distribution from the Ising model.
        Q_unnormalized = binom.*exp(states*param_list);

        % Normalize this Q.
        % This gives us a probability distribution Q(k), the probability of
        % finding a population with count k.
        Q = Q_unnormalized/sum(Q_unnormalized);

        % Find non-zero entries as we will be taking the log.
        ind_p = find(P ~= 0);
        ind_q = find(Q ~= 0);

        % Calculate the DKL.
        DKL(k) = -Q(ind_q)'*log2(Q(ind_q)) + P(ind_p)'*log2(P(ind_p));
    end

    % Sort rho to get it in increasing order.
    [rho_sorted,IX] = sort(rho);
    % Match corresponding DKLs.
    DKL_sorted = DKL(IX);

    DKL_averaged = zeros(length(M{j})/avnum,1);
    rho_averaged = zeros(length(M{j})/avnum,1);

    % Average over avnum values to get a smoother plot.
    for i=1:length(M{j})/avnum
        DKL_averaged(i) = mean(DKL_sorted((avnum*(i-1)+1):avnum*i));
        rho_averaged(i) = mean(rho_sorted((avnum*(i-1)+1):avnum*i));
    end

    % Sort once more.
    [rho_final,IX] = sort(rho_averaged);
    DKL_final = DKL_averaged(IX);

    plot(rho_final,DKL_final,'color',CM(j,:))
    hold on
end

xlabel('Correlation coefficient \rho','fontsize',16)
ylabel('KL-divergence \Delta_h','fontsize',16)

axis([0 0.3 0 10e-3])

toc