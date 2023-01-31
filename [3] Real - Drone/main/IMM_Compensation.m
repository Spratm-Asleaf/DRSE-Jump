%{
Online supplementary materials of the paper titled 
"Distributionally Robust State Estimation for Jump Linear Systems"
Authored By: Shixiong Wang
From the Institute of Data Science, National University of Singapore

@Author: Shixiong Wang
@Date: 16 May 2022
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE-Jump
%}

function [x_hat, mu_hat] = IMM_Compensation(F, G, H, Q, R, Pi, Y, a, x_hat0, P0, mu_hat0)
%{
    This is an implementation of 
    Zhao (2016).
    Recursive estimation for Markov jump linear systems with unknown transition probabilities: A compensation approach. 
    Journal of Franklin Institute.
%}

n = length(x_hat0);
[~, episode_length] = size(Y);
[N,~] = size(Pi);

% merged posterior estimate
x_hat = zeros(n, episode_length);

% individual estiamte; we have "N" models
X = zeros(n, N);
P = cell(1, N);
for j = 1:N         % init for every candidate model
    X(:, j) = x_hat0;
    P{j} = P0;
end

X_temp = zeros(n, N);
P_temp = cell(1, N);
    
% estimated model probability
mu_hat = zeros(N, episode_length);

% filter starts
for k = 1:episode_length
    for j = 1:N
        % Step 1.1
        if k == 1
            mu_ij = Pi(:,j) .* mu_hat0;
        else 
            mu_ij = Pi(:,j) .* mu_hat(:, k-1);
        end
        sum_mu_ij = sum(mu_ij);
        mu_ij = mu_ij/sum(mu_ij);
        
        % Step 1.2
        X_temp(:, j) = zeros(n, 1);
        P_temp{j} = zeros(size(P0));
        for i = 1:N
            X_temp(:, j) = X_temp(:, j) + mu_ij(i) * X(:, i);
        end
        for i = 1:N             % do not merge this two "for" loops
            P_temp{j} = P_temp{j} + mu_ij(i) * (P{i} + (X(:, i) - X_temp(:, j))*(X(:, i) - X_temp(:, j))');
        end

        % Step 1.3
        X_temp(:, j) = F*X_temp(:, j) + G*a(:, j);
        P_temp{j} = F*P_temp{j}*F' + G*Q*G';
        
        % Step 1.4
        r = Y(k) - H*X_temp(:, j);
        S = H*P_temp{j}*H' + R;
        K = P_temp{j}*H'*S^-1;
        X_temp(:, j) = X_temp(:, j) + K*r;
        P_temp{j} = P_temp{j} - P_temp{j}*H'*S^(-1)*H*P_temp{j};
        
        % Step 1.5 ~ Step 1.7
        mu_hat(j, k) = sum_mu_ij * normpdf(r, 0, S); % N.B.: This "normpdf" function is different from the MATLAB's built-in version
    end
    
    % update X and P values
    X = X_temp;
    P = P_temp;
    
    % Step 1.7, normalization of mu_hat
    mu_hat_sum = sum(mu_hat(:, k));
    for j = 1:N      % do not merge this two "for" loops
        mu_hat(j, k) = mu_hat(j, k)/mu_hat_sum;
    end

    %%{{ BEGIN of Zhao's compensation; See the note at the end of this file.
        % Eq. (30) therein
        chi = max(mu_hat(:, k));
        kappa = min(mu_hat(:, k));
        sigma = chi + kappa;

        if N > 2
            if sigma > 2/N
                eta = (chi - 1/N)/abs(kappa - 1/N);
            else
                eta = 1;
            end
        elseif N == 2
            alpha = 0.6;             % as suggested in "Example 1" therein
            if chi - alpha >= 0.1
                eta = 10*(chi - alpha);
            else
                eta = 1;
            end
        end
        
        % Eq. (28) therein
        mu_hat(:, k) = (mu_hat(:, k)).^eta;
        mu_hat_sum = sum(mu_hat(:, k));
        mu_hat(:, k) = mu_hat(:, k)/mu_hat_sum;
    %%}} End of Zhao's compensation

    % Step 2
    x_hat(:, k) = zeros(n, 1);
    for j = 1:N
        x_hat(:, k) = x_hat(:, k) + mu_hat(j, k) * X(:, j);
    end
    
    % No need to calculate mixted posterior estimation error covaraince
%     P_temp = zeros(size(P0));
%     for j = 1:N             % do not merge this two "for" loops
%         P_temp = P_temp + mu_hat(j, k) * (P{j} + (X(:, j) - x_hat(:, k))*(X(:, j) - x_hat(:, k))');
%     end
end

%%%%%%%%%%%%%%  Note %%%%%%%%%%%%%%%%%%
%{
% "Example 1" inside of (Zhao 1016; Journal of Franklin Institute)
% Readers may use the codes below to verify Zhao's claims in his "Example 1" 
% input = [0.7, 0.1, 0.1, 0.1]
% output = [0.991, 0.003, 0.003, 0.003, 0.003]
% The codes gives the same answer as claimed by Zhao in (Zhao 2016; Journal of Franklin Institute)
% Hence, my codes to implement Zhao's compensation method is loyal to Zhao.

mu_hat = [0.7 0.1 0.1 0.1];
N = 4;
chi = max(mu_hat);
kappa = min(mu_hat);
sigma = chi + kappa;

if N > 2
    if sigma > 2/N
        eta = (chi - 1/N)/abs(kappa - 1/N);
    else
        eta = 1;
    end
elseif N == 2
    alpha = 0.6;             % as suggested in "Example 1" inside
    if chi - alpha >= 0.1
        eta = 10*(chi - alpha);
    else
        eta = 1;
    end
end

mu_hat = mu_hat.^eta;
mu_hat_sum = sum(mu_hat);
mu_hat = mu_hat/mu_hat_sum;
%}