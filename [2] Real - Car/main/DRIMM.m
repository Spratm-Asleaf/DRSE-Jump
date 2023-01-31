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

function [x_hat, mu_hat] = DRIMM(Algo, F, G, H, Q, R, Pi, Y, a, x_hat0, P0, mu_hat0)
    n = length(x_hat0);
    [~, episode_length] = size(Y);
    [N,~] = size(Pi);

    % merged posterior estimate
    x_hat = zeros(n, episode_length);

    % individual estiamte; we have "N" models
    X = zeros(n, N);
    P = cell(1, N);
    for j = 1:N         % init for every candidate model
        P{j} = P0;% * randi(5);
        X(:, j) = x_hat0;% + chol(P{j}) * randn(n, 1);
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

        % update X and P values
        X = X_temp;
        P = P_temp;

        % Step 1.7, normalization of mu_hat
        mu_hat_sum = sum(mu_hat(:, k));
        for j = 1:N      % do not merge this two "for" loops
            mu_hat(j, k) = mu_hat(j, k)/mu_hat_sum;
        end

        % Step 2, Distributionally Robust Start
        if true
            Algo = strrep(Algo,' ','');
            Algo = lower(Algo);
            muBar = mu_hat(:, k);
            switch Algo
                % Robustify over "\mu"
                case 'mu'
                    mu_hat(:, k) = RobustifyOverMu(muBar, X, P);
                        % There is another implementation:
                        % RobustifyOverMu2(muBar, X, P)
                        % But the second version is numerically not stalbe.
                        % See details inside the function "RobustifyOverMu2"
                % Robustify over "P" 
                case 'p'
                    P = RobustifyOverP(P);
                % Robustify over both "\mu" and "P"
                case 'mu+p'
                    % Do note that modify "P" first, modify "mu" next.
                    % Cannot reverse this order.
                    P = RobustifyOverP(P);
                    mu_hat(:, k) = RobustifyOverMu(muBar, X, P);
                otherwise
                    error('DRIMM:: Specified Algorithm Not Found! Check Spelling.');
            end
        end

        % Step 3
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
    
function mu = RobustifyOverMu(muBar, X, P)
% This function implements Algorithm 3 in the online supplementary materials.
%
% We are required to solve a nonlinear root-finding problem.
% Readers can solve the associated nonlinear root-finding problem.
%
% However, below is another trick to walk it around.
% The basic idea is to solve (30) but the objective is replaced with (23)
% We iterate:
% (1) Fix "mu1" to find weighted mean, and maximize over "mu", i.e.,
%           maximize sum_j {mu_j * [P_j + (x_j - mu1*x)(x_j - mu1*x)']}.
%     Suppose "mu2" solves the above problem.
% (2) Replace mu1 with mu2, and then go to step (1), until the iteration process converges.
% Note that (1) has a linear objective and hence easy to be solved with KKT conditions (i.e., Lagrangian duality method).
% I did not discuss this trick in the main body of the paper because I cannot rigorously prove its convergence, although it empirically works well.
% I am happy to see if anyone can prove this method.

    outerLoop = 100;
    theta_0 = 0.01;
    
    mu = muBar;
    
    N = length(muBar);
    
    isFirstStepIn = true;
    order = 0;

    %% Iteration starts
    lastCost = 0;
    while outerLoop >= 0
        outerLoop = outerLoop - 1;

        %% Step (2): Find weightd mean
        X_Avg = X*mu;
        c = zeros(N, 1);
        for j = 1:N
            c(j) = trace(P{j} + (X(:, j) - X_Avg)*(X(:, j) - X_Avg)');
        end
        
        % N.B.: Cannot execute this statement in every WHILE loop
        % Because the scale of "c" can be changed in the successive loops
        if isFirstStepIn
            % Re-scale
            % The numerical scale of A and b is very large and therefore re-scaling is important.
            order = 0;
            max_value = max(c);
            while true
                if max_value * 10^(-order) > 1
                    order = order + 1;
                else
                    break;
                end
            end
            
            isFirstStepIn = false;      % N.B., Very Important. Do not Delete.
        end
        
        c = c*10^(-order);

        currentCost = c'*mu;
        if abs(currentCost - lastCost) < 1e-3
            break;
        end
        lastCost = currentCost;

        %% Step (1): Find maximizer
        % Maximize over "mu" using Lagrangian duality (i.e., KKT conditions)
        innerLoop = 200;            % Do not set this value to be extremely large to pursue the optimality; 200 is good for this experiment.
        lambda0 = 5;                % Any positive value is ok to initialize "lambda0"
        lambda1 = -lambda0;         % Any real value is ok to initialize "lambda1"
        alpha = 0.1;               % Step size

        p = zeros(N, 1);            % Improved value of "mu" is "p"
        while innerLoop >= 0
            innerLoop = innerLoop - 1;
            for j = 1:N
                % Use KKT conditions for the maximization problem in Step (1) to see the following formula.
                p(j) = exp((c(j) - lambda1 - lambda0 + lambda0*log(muBar(j)))/lambda0);
            end

            grad_lambda0 = theta_0 - (p'*log(p) - p'*log(muBar));
            grad_lambda1 = 1 - sum(p);

            lambda0 = lambda0 - alpha * grad_lambda0;
            if lambda0 <= 0
                lambda0 = 0;
            end
            lambda1 = lambda1 - alpha * grad_lambda1;

            if abs(sum(p) - 1) < 1e-2
                break;
            end
        end

        mu = p/sum(p);
    end
end

function mu = RobustifyOverMu2(muBar, X, P)
% This function implements Algorithm 3 in the online supplementary materials
%
% This algorithm is numerically not stable.
% When the augmented matrix "[A b]" is singular, then there is no improvement of the objective over "mu".
% Whether this condition will be triggled or not dependes on the characteristics of problems.
% i.e., finally, "mu = muBar".
    N = length(muBar);
    A = zeros(N, N);
    b = zeros(N, 1);
    for i = 1:N
        for j = 1:N
            A(i, j) = (X(:, i))' * X(:, j);
        end
        b(i) = trace(P{i}) + (X(:, i))' * X(:, i);
    end

    % Re-scale
    % The numerical scale of A and b is very large and therefore re-scaling
    % is important.
    order = 0;
    max_value = max(b);
    while true
        if max_value * 10^(-order) > 1
            order = order + 1;
        else
            break;
        end
    end

    A = A * 10^(-order);
    b = b * 10^(-order);

    theta_0 = 0.1;

    [mu,~,~] = fmincon(@(mu)mu'*A*mu - b'*mu,...
        muBar,...           % init start point
        [], [],...          % linear inequality 
        ones(1, N), 1,...   % linear equality
        zeros(N, 1),...     % lower bound
        ones(N, 1),...      % upper bound
        @(mu)NonlinearConstraint(mu, muBar, theta_0)...
    );
end

function [c, ceq] = NonlinearConstraint(mu, muBar, theta0)
    c = (mu - muBar)'*(mu - muBar) - theta0;    % Norm Constraint
%     c = mu.*(ln(mu) - ln(muBar))  - theta0;     % KL Divergence
    ceq = [];
end

function P = RobustifyOverP(P)
    N = length(P);
    
    theta_j = 1.25;
    
    for j = 1:N
        if j == 1
            continue; % the first model is exact: the value of "a(1)" is zero for both true and uncertain models.
        end
        P{j} = P{j} * theta_j;
    end
end