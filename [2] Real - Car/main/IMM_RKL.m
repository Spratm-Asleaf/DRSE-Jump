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

function [x_hat, mu_hat] = IMM_RKL(F, G, H, Q, R, Pi, Y, a, x_hat0, P0, mu_hat0)
%{
    This is an implementation of 
    Orguner (2006).
    An online sequential algorithm for the estimation of transition probabilities for jump Markov linear systems.
    Automatica.
%}

% Orguner (2006) claimed to use equal initial weights
% mu_hat0 = ones(3,1)/3;

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

%%{{ BEGIN of Orguner's recursive KL
    % estimated model transition probability matrix
    % Eq. (80) therein
    Pi_hat = [
        0.33 0.33 0.34
        0.33 0.33 0.34
        0.33 0.33 0.34
    ];                  % I found this value useless
    Pi_hat = Pi;        % So I tried the nominal value
%%}} End of Orguner's recursive KL

% filter starts
for k = 1:episode_length
    for j = 1:N
        % Step 1.1
        if k == 1
            mu_ij = Pi_hat(:,j) .* mu_hat0;
        else 
            mu_ij = Pi_hat(:,j) .* mu_hat(:, k-1);
        end
        sum_mu_ij = sum(mu_ij);
        
        if sum_mu_ij < 1e-6
            mu_ij = mu_ij/1e-6;
        else
            mu_ij = mu_ij/sum(mu_ij);
        end
        
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
        X_temp(:, j) = F*X_temp(:, j) + G*a(j);
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
    
    %%{{ BEGIN of Orguner's recursive KL
        if k >= 2
            % calculate the denominator in the TPM's update rule
            den = 0;
            for i = 1:N
                for j = 1:N
                    % Eq. (72) therein
                    y_ij = H*(F*X(:,i) + G*a(j));
                    % EQ. (73) therein
                    S_ij = H*(F*P{j}*F'+G*Q*G')*H' + R;
                    % calculate den
                    den = den + normpdf(Y(k), y_ij, S_ij)*Pi_hat(i,j)*mu_hat(i, k-1);
                end
            end
            
            if den < 1e-6
                den = 1e-6;
            end

            % Update TPM, i.e., Ph_hat
            epsilon = 0.01; % used in its experiments
            for i = 1:N
                for j = 1:N
                    % Eq. (72) therein
                    y_ij = H*(F*X(:,i) + G*a(j));
                    % Eq. (73) therein
                    S_ij = H*(F*P{j}*F'+G*Q*G')*H' + R;
                    % Update rule below Eq. (73)
                    Pi_hat(i,j) = Pi_hat(i,j) + epsilon*(normpdf(Y(k), y_ij, S_ij)*mu_hat(i, k-1))/den;
                end
            end

            % Projection
            for i = 1:N
                Pi_hat(i, :) = Project(Pi_hat(i, :), N);
            end
        end
    %%}} End of Orguner's recursive KL
    
    % Step 1.7, normalization of mu_hat
    mu_hat_sum = sum(mu_hat(:, k));
    for j = 1:N      % do not merge this two "for" loops
        mu_hat(j, k) = mu_hat(j, k)/mu_hat_sum;
    end
    
    % update X and P values
    X = X_temp;
    P = P_temp;
    
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
end


%% Project a vector into a probability simplex
% The following function is verified to be correct because its result is
% consistent with the projection method in 
% https://www.mathworks.com/matlabcentral/fileexchange/30332-projection-onto-simplex?s_tid=mwa_osa_a
function p = Project(x, N)
% This method was provided in Appendix A of Orguner (2006)
    if N == 1
        p = min(x, 1);
        return;
    end
    
    p = zeros(1, N);
    e = (sum(x) - 1)/N;
    [min_val, jj] = min(x);
    if e <= min_val
        p = x - e*ones(1,N);
    else
        y = [x(1:jj-1) x(jj+1:N)];
        m = Project(y, N-1);
        if 1 <= jj-1
            p(1:jj-1) = m(1:jj-1);
        else
            p(1:jj-1) = [];
        end
        p(jj) = 0;
        
        if jj <= N-1
            p(jj+1:N) = m(jj:N-1);
        else
            p(jj+1:N) = [];
        end
    end
end
