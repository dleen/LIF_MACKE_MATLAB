clear all;

% Choose between Matlab optimization = 1
% or minimize.m minimization function = 2.
alg_choice = 1;

% Number of neurons.
N = 100;

% Load data from a file. 
% The file.dat has the form:
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
% We keep mu and rho fixed over trials. We have multiple
% trials so that we can get an average and smooth out
% some of the inherent noise.
M = importdata('figure_2a_0.64.dat',' ');

% PP is a matrix where each row is a prob dist.
PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

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

% Generate a list of binomial coefficients.
binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

% An initial guess at the parameters to be minimized:
param_list_init = [0;0];


% Take the mean of the input p(x).
% i.e we average over the trials
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 1
% p(x_1) p(x_2) ... p(x_N) mu rho % trial number 2
P = mean(PP,1)';

% Output mean firing rate and corr coefficient.
mean_mu = mean(mu)
mean_rho = mean(rho)

% Calculate the means and moments i.e. the quantities
% sum(k*P_k) and sum(k^2*P_k).
mean_feature = P'*states;

% Choose the algorithm to use.
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

% Use the parameters found via minimization to calculate the probability
% distribution from the Ising model.
% The binomial prefactor is explained in Macke 2011. 
Q_unnormalized = binom.*exp(states*param_list);

% Normalize this Q.
% This gives us a probability distribution Q(k), the probability of
% finding a population with count k.
Q = Q_unnormalized/sum(Q_unnormalized);

% Plot output.
semilogy(Q,'b');
hold on
semilogy(P,'r');
axis([0 N 1e-4 1]);
legend('Ising Model Prediction','Numerical Simulation')