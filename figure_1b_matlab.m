clear all;

alg_choice = 1;

N = 5;

tic
M = importdata('figure_1b_0.65.dat',' ');

PP = M(:,1:N+1);
mu = M(:,N+2);
rho = M(:,N+3);

states = [(0:N)',(0:N)'.^2];

binom = zeros(N+1,1);
for i=0:N
    binom(i+1) = nchoosek(N,i);
end

% An initial guess:
param_list_init = [0;0];

options = optimset('GradObj','on','LargeScale','on',...
        'Display','off',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-7,'TolX',1e-7);

for k=1:length(M)
    P = (PP(k,:))';

    % Calculate the means and moments i.e. the quantities
    % sum(k*P_k) and sum(k^2*P_k).
    mean_feature = P'*states;

    % Eric's minimization function, based on conjugate gradient.
     if alg_choice == 1
    % MATLAB Optimization Toolbox minimization function.
    % Uses Optimization toolbox unconstrained minimization function.
        [param_list] =...
            fminunc(@(x)neg_log_like_binom...
            (x,states,mean_feature,P,binom),param_list_init,options);
     end

    % Use the parameters found via minimization to calculate the probability
    % distribution from the Ising model.
    Q_unnormalized = binom.*exp(states*param_list);

    % Normalize this Q.
    % This gives us a probability distribution Q(k), the probability of
    % finding a population with count k.
    %QQ(k,:) = Q_unnormalized/sum(Q_unnormalized);
    Q = Q_unnormalized/sum(Q_unnormalized);

    
%     ind_p = find(P > 1e-15);
%     ind_q = find(Q > 1e-15);
    
     ind_p = find(P ~= 0);
     ind_q = find(Q ~= 0);

    DKL = -Q(ind_q)'*log2(Q(ind_q)) + P(ind_p)'*log2(P(ind_p));
    DKLlist(k) = DKL;
end

DKLlist = DKLlist';

avnum = 50;

for i=1:length(M)/avnum
    DKLlistnew(i) = mean(DKLlist((avnum*(i-1)+1):avnum*i));
    rhonew(i) = mean(rho((avnum*(i-1)+1):avnum*i));
end
[B,IX] = sort(rhonew);
rho = B;
DKLlist = DKLlistnew(IX);

% [B,IX] = sort(rho);
% rho = B;
% DKLlist = DKLlist(IX);
toc

hold on
plot(rho,DKLlist,'r')
axis([0 0.3 0 15e-3])