clear all;

alg_choice = 1;

N = 100;

tic
M = importdata('figure_2a_0.7.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);


states = [(0:N)',(0:N)'.^2];

%QQ = zeros(size(PP));

binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

% An initial guess:
param_list_init = [0;0];

%for k=1:size(M,1)
%    P = (PP(k,:))';

P = mean(PP)';

    % Calculate the means and moments i.e. the quantities
    % sum(k*P_k) and sum(k^2*P_k).
    mean_feature = P'*states;

    if alg_choice == 1
    % MATLAB Optimization Toolbox minimization function.
    % Uses Optimization toolbox unconstrained minimization function.
    options = optimset('GradObj','on','LargeScale','on',...
        'Display','final-detailed',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
        [param_list,fval,output] =...
            fminunc(@(x)neg_log_like_binom...
            (x,states,mean_feature,P,binom),param_list_init,options);
elseif alg_choice == 2
    % Eric's minimization function, based on conjugate gradient.
    % Previous value '200' instead of '400'
    param_list = minimize(param_list_init,...
          'neg_log_like_binom',400,...
          states,mean_feature,P,binom);
end
    

    % Use the parameters found via minimization to calculate the probability
    % distribution from the Ising model.
    Q_unnormalized = binom.*exp(states*param_list);

    % Normalize this Q.
    % This gives us a probability distribution Q(k), the probability of
    % finding a population with count k.
    %QQ(k,:) = Q_unnormalized/sum(Q_unnormalized);
    Q = Q_unnormalized/sum(Q_unnormalized);

 %   QQ(k,:) = Q;
    
%end

semilogy(Q,'b');
hold on
%semilogy(P,'b');
axis([0 N 1e-4 1]);