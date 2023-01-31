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

clc;
clear all;
close all;

global scenario_index;

% If you want to exactly recover the results reported in my paper, keep the following statement.
% You can also try other values to see what will happen, if interested.
if true
    rng(46);     % Just because I like the number 46: My mom is 46 years old this year.
end

NumMonteCarlo = 50;                    % how many independent episodes to run
EpisodeLength = 1000;                   % how many time steps included in each episode
% (KF) - Kalman Filter with known true system dynamics
g_err_KF = cell(1, NumMonteCarlo);      % "g" for "global"
g_time_KF = zeros(1, NumMonteCarlo);
g_RMSE_KF = zeros(1, NumMonteCarlo);
% (IMM-T) - IMM with true model transition probability matrix
g_err_IMM_true = cell(1, NumMonteCarlo);
g_time_IMM_true = zeros(1, NumMonteCarlo);
g_RMSE_IMM_true = zeros(1, NumMonteCarlo);
% (IMM-N) - IMM with nominal (i.e., user-designed) model transition probability matrix
g_err_IMM_nominal = cell(1, NumMonteCarlo);
g_time_IMM_nominal = zeros(1, NumMonteCarlo);
g_RMSE_IMM_nominal = zeros(1, NumMonteCarlo);
% (IMM-C) - IMM with compensation; Zhao (2016) "Recursive estimation for Markov jump linear systems with unknown transition probabilities: A compensation approach"
g_err_IMM_compensation = cell(1, NumMonteCarlo);
g_time_IMM_compensation = zeros(1, NumMonteCarlo);
g_RMSE_IMM_compensation = zeros(1, NumMonteCarlo);
% (IMM-B) - IMM with Bayesian estimation; Jilkov (2004) "Online Bayesian estimation of transition probabilities for Markovian jump systems"
g_err_IMM_Bayesian = cell(1, NumMonteCarlo);
g_time_IMM_Bayesian = zeros(1, NumMonteCarlo);
g_RMSE_IMM_Bayesian = zeros(1, NumMonteCarlo);
% (IMM-M) - IMM with maximum likelihood (specifically, Recursive Kullback–Leibler); Orguner (2006) "An online sequential algorithm for the estimation of transition probabilities for jump Markov linear systems"
g_err_IMM_MaximumLikelihood = cell(1, NumMonteCarlo);
g_time_IMM_MaximumLikelihood = zeros(1, NumMonteCarlo);
g_RMSE_IMM_MaximumLikelihood = zeros(1, NumMonteCarlo);
% (IMM-RS) - risk-sensitive IMM; Orguner (2008) "Risk-sensitive filtering for jump Markov linear systems"
g_err_IMM_RiskSensitive = cell(1, NumMonteCarlo);
g_time_IMM_RiskSensitive = zeros(1, NumMonteCarlo);
g_RMSE_IMM_RiskSensitive = zeros(1, NumMonteCarlo);
% (IMM-R) - IMM using the distributionally robust Kalman filter; Wang (2021) "Robust state estimation for linear systems under distributional uncertainty"
g_err_IMM_RobustKF = cell(1, NumMonteCarlo);
g_time_IMM_RobustKF = zeros(1, NumMonteCarlo);
g_RMSE_IMM_RobustKF = zeros(1, NumMonteCarlo);
% (DRIMM) - Distributionally Robust IMM; This paper
g_err_DRIMM = cell(1, NumMonteCarlo);
g_time_DRIMM = zeros(1, NumMonteCarlo);
g_RMSE_DRIMM = zeros(1, NumMonteCarlo);

%% main filtering algorithms
for episode = 1:NumMonteCarlo
    clc;
    disp(['Episode: ' num2str(episode) '/' num2str(NumMonteCarlo)]);
    
    %% Load simulation scenario (i.e., the true system data)
    % a : nominal (i.e., user-designed) value of "sys.a", which might be different from the true value "sys.a"
    % Pi: nominal (i.e., user-designed) value of "sys.Pi", which might be different from the true value "sys.Pi"
    scenario_index = 6;
    switch scenario_index   
        case 1      % Table I
            true_Pi = [
                0.8     0.1     0.1
                0.1     0.8     0.1
                0.1     0.1     0.8
            ];
        
            sys = LoadScenario1(EpisodeLength, true_Pi);
            
            a = sys.a;  
            
            Pi = [
                0.9     0.05    0.05
                0.05    0.9     0.05
                0.05    0.05    0.9
            ];
        
        case 2      % Table II 
            true_Pi = [
                0.1     0.1     0.8
                0.1     0.8     0.1
                0.1     0.1     0.8
            ];
        
            sys = LoadScenario1(EpisodeLength, true_Pi);
            
            a = sys.a;  
            
            Pi = [
                0.8     0.1     0.1
                0.1     0.8     0.1
                0.1     0.1     0.8
            ];
            
        case 3      % Table III
            sys = LoadScenario2(EpisodeLength);

            a = sys.a(1:3);             % nominal "a" contains the first three values of true "sys.a"
               
            Pi = [                      % nominal "Pi" is 3*3 but true "sys.Pi" is 4*4
                0.6     0.2     0.2
                0.2     0.6     0.2
                0.2     0.2     0.6
            ];

        case 4      % Table IV
            true_Pi = [
                0.8     0.1     0.1
                0.1     0.8     0.1
                0.1     0.1     0.8
            ]; 
        
            sys = LoadScenario1(EpisodeLength, true_Pi);
            
            a = 0.5*sys.a;
            
            Pi = sys.Pi;
            
        case 5      % Table V
            true_Pi = [
                0.8     0.1     0.1
                0.1     0.8     0.1
                0.1     0.1     0.8
            ]; 
        
            sys = LoadScenario1(EpisodeLength, true_Pi);
            
            a = 0.5*sys.a;  
            
            Pi = [
                0.6     0.2     0.2
                0.2     0.6     0.2
                0.2     0.2     0.6
            ];
        
        case 6      % Table VI
            true_Pi = [
                0.8     0.1     0.1
                0.1     0.8     0.1
                0.1     0.1     0.8
            ]; 
        
            sys = LoadScenario1(EpisodeLength, true_Pi);
            
            a = sys.a;  
            
            Pi = sys.Pi;
    end
    
    %% Initialization
    x_hat0  = sys.init.x0;
    P0      = sys.init.P0;
    if scenario_index == 3
        mu_hat0 = [0.8;0.1;0.1];    % nominal "mu_hat0" is 3*1 but true "mu_hat0" is 4*1
    else
        mu_hat0 = sys.init.mu0;
    end

    %% In which dimension do we calculate error statistics?
    err_dim = 1:2;    % 1 (position), 2 (speed), 1:2 (all)

    %% KF (Kalman filter using the true model trajectory; This is the optimal method)
    % KF
    % N.B.: The true model trajectory is, in princple and in practice, unknown. But in simulations, we know it exactly.
    tic
    [x_hat_KF, ~] = KF(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.y, sys.j, sys.a, x_hat0, P0);
    time_KF = toc/sys.K;
    %disp(['Avg Time of KF: ' num2str(time_KF)]);
    err_KF = sum((sys.x(err_dim, :) - x_hat_KF(err_dim, :)).^2,1);
    RMSE_KF = sqrt(mean(err_KF));
    %disp(['---------------------RMSE of KF: ' num2str(RMSE_KF)]);
    
    % save data for each episode
    g_err_KF{episode} = err_KF;
    g_time_KF(episode) = time_KF;
    g_RMSE_KF(episode) = RMSE_KF;

    %% Standard IMM (True Pi & True a)
    % IMM-T
    % N.B.: We know the true TPM, i.e., Pi, in simulations.
    tic
    [x_hat_IMM_true, mu_hat_IMM_true] = IMM(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, sys.init.mu0);
    time_IMM_true = toc/sys.K;
    %disp(['Avg Time of IMM (True): ' num2str(time_IMM_true)]);
    err_IMM_true = sum((sys.x(err_dim, :) - x_hat_IMM_true(err_dim, :)).^2,1);
    RMSE_IMM_true = sqrt(mean(err_IMM_true));
    %disp(['---------------------RMSE of IMM (True): ' num2str(RMSE_IMM_true)]);
    
    % save data for each episode
    g_err_IMM_true{episode} = err_IMM_true;
    g_time_IMM_true(episode) = time_IMM_true;
    g_RMSE_IMM_true(episode) = RMSE_IMM_true;
    
    %% Standard IMM (Nominal Pi and Nominal a)
    % IMM-N
    tic
    [x_hat_IMM_nominal, mu_hat_IMM_nominal] = IMM(sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_nominal = toc/sys.K;
    %disp(['Avg Time of IMM (Nominal): ' num2str(time_IMM_nominal)]);
    err_IMM_nominal = sum((sys.x(err_dim, :) - x_hat_IMM_nominal(err_dim, :)).^2,1);
    RMSE_IMM_nominal = sqrt(mean(err_IMM_nominal));
    %disp(['---------------------RMSE of IMM (Nominal): ' num2str(RMSE_IMM_nominal)]);
    
    % save data for each episode
    g_err_IMM_nominal{episode} = err_IMM_nominal;
    g_time_IMM_nominal(episode) = time_IMM_nominal;
    g_RMSE_IMM_nominal(episode) = RMSE_IMM_nominal;
    
    %% Compensation IMM; Zhao Shunyi (2016) Journal of Franklin Institute
    % IMM-C
    % Recursive estimation for Markov jump linear systems with unknown transition probabilities: A compensation approach. 
    tic
    [x_hat_IMM_compensation, mu_hat_IMM_compensation] = IMM_Compensation(sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_compensation = toc/sys.K;
    %disp(['Avg Time of IMM (Compensation): ' num2str(time_IMM_compensation)]);
    err_IMM_compensation = sum((sys.x(err_dim, :) - x_hat_IMM_compensation(err_dim, :)).^2,1);
    RMSE_IMM_compensation = sqrt(mean(err_IMM_compensation));
    %disp(['---------------------RMSE of IMM (Compensation): ' num2str(RMSE_IMM_compensation)]);
    
    % save data for each episode
    g_err_IMM_compensation{episode} = err_IMM_compensation;
    g_time_IMM_compensation(episode) = time_IMM_compensation;
    g_RMSE_IMM_compensation(episode) = RMSE_IMM_compensation;
	
	
    %% IMM with Bayesian estimation; Jilkov and Li (2004) IEEE Trans on Signal Processing
    % IMM-B
    % Online Bayesian estimation of transition probabilities for Markovian jump systems
    tic
    [x_hat_IMM_Bayesian, mu_hat_IMM_Bayesian] = IMM_Bayesian(sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_Bayesian = toc/sys.K;
    %disp(['Avg Time of IMM (Bayesian): ' num2str(time_IMM_Bayesian)]);
    err_IMM_Bayesian = sum((sys.x(err_dim, :) - x_hat_IMM_Bayesian(err_dim, :)).^2,1);
    RMSE_IMM_Bayesian = sqrt(mean(err_IMM_Bayesian));
    %disp(['---------------------RMSE of IMM (Bayesian): ' num2str(RMSE_IMM_Bayesian)]);
    
    % save data for each episode
    g_err_IMM_Bayesian{episode} = err_IMM_Bayesian;
    g_time_IMM_Bayesian(episode) = time_IMM_Bayesian;
    g_RMSE_IMM_Bayesian(episode) = RMSE_IMM_Bayesian;
    
    %% IMM with Maximum Likelihood estimation; The Recursive-Kullback--Leibler (RKL) method; Orguner (2006)
    % IMM-M
    tic
    [x_hat_IMM_MaximumLikelihood, mu_hat_IMM_MaximumLikelihood] = IMM_RKL(sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_MaximumLikelihood = toc/sys.K;
    %disp(['Avg Time of IMM (Maximum Likelihood, RKL): ' num2str(time_IMM_RKL)]);
    err_IMM_MaximumLikelihood = sum((sys.x(err_dim, :) - x_hat_IMM_MaximumLikelihood(err_dim, :)).^2,1);
    RMSE_IMM_MaximumLikelihood = sqrt(mean(err_IMM_MaximumLikelihood));
    %disp(['---------------------RMSE of IMM (RKL): ' num2str(RMSE_IMM_RKL)]);
    
    % save data for each episode
    g_err_IMM_MaximumLikelihood{episode} = err_IMM_MaximumLikelihood;
    g_time_IMM_MaximumLikelihood(episode) = time_IMM_MaximumLikelihood;
    g_RMSE_IMM_MaximumLikelihood(episode) = RMSE_IMM_MaximumLikelihood;
    
    %% Risk-Sensitive IMM; Orguner (2008)
    % IMM-RS
    tic
    [x_hat_IMM_RiskSensitive, mu_hat_IMM_RiskSensitive] = IMM_RS(sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_RiskSensitive = toc/sys.K;
    %disp(['Avg Time of IMM (Risk-Sensitive): ' num2str(time_IMM_RS)]);
    err_IMM_RiskSensitive = sum((sys.x(err_dim, :) - x_hat_IMM_RiskSensitive(err_dim, :)).^2,1);
    RMSE_IMM_RiskSensitive = sqrt(mean(err_IMM_RiskSensitive));
    %disp(['---------------------RMSE of IMM (RS): ' num2str(RMSE_IMM_RS)]);
    
    % save data for each episode
    g_err_IMM_RiskSensitive{episode} = err_IMM_RiskSensitive;
    g_time_IMM_RiskSensitive(episode) = time_IMM_RiskSensitive;
    g_RMSE_IMM_RiskSensitive(episode) = RMSE_IMM_RiskSensitive;
    
    %% IMM using Robust Kalman Filter; Wang (2021)
    % IMM-R
    tic
    [x_hat_IMM_RobustKF, mu_hat_IMM_RobustKF] = DRIMM('P', sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_IMM_RobustKF = toc/sys.K;
    %disp(['Avg Time of IMM (R): ' num2str(time_IMM_R)]);
    err_IMM_RobustKF = sum((sys.x(err_dim, :) - x_hat_IMM_RobustKF(err_dim, :)).^2,1);
    RMSE_IMM_RobustKF = sqrt(mean(err_IMM_RobustKF));
    %disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_IMM_R)]);
    
    % save data for each episode
    g_err_IMM_RobustKF{episode} = err_IMM_RobustKF;
    g_time_IMM_RobustKF(episode) = time_IMM_RobustKF;
    g_RMSE_IMM_RobustKF(episode) = RMSE_IMM_RobustKF;
    
    %% Distributionally Robust IMM; This paper
    % DRIMM
    tic
    % Possible Values of "Algo": 'mu', 'P', 'mu+P'
    if scenario_index <= 2      % only TPM is uncertain and therefore model weights are uncertain
        Algo = 'mu'; 
    elseif scenario_index == 3  % model set not complete: using three among five
        Algo = 'mu+P';
    elseif scenario_index == 4  % only model uncertainty 
        Algo = 'mu+P';
    elseif scenario_index == 5  % model uncertainty + uncertain TPM
        Algo = 'mu+P';
    elseif scenario_index == 6  % no uncertainty 
        Algo = 'mu+P'; 
    end
    [x_hat_DRIMM, mu_hat_DRIMM] = DRIMM(Algo, sys.F, sys.G, sys.H, sys.Q, sys.R, Pi, sys.y, a, x_hat0, P0, mu_hat0);
    time_DRIMM = toc/sys.K;
    %disp(['Avg Time of DRIMM: ' num2str(time_DRIMM)]);
    err_DRIMM = sum((sys.x(err_dim, :) - x_hat_DRIMM(err_dim, :)).^2,1);
    RMSE_DRIMM = sqrt(mean(err_DRIMM));
    %disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_DRIMM)]);
    
    % save data for each episode
    g_err_DRIMM{episode} = err_DRIMM;
    g_time_DRIMM(episode) = time_DRIMM;
    g_RMSE_DRIMM(episode) = RMSE_DRIMM;
end

%% Scenario
if scenario_index == 1      % only uncertain TPM
    prompt = 'Only uncertainty in TPM - Case 1';
elseif scenario_index == 2  % only uncertain TPM
    prompt = 'Only uncertainty in TPM - Case 2';
elseif scenario_index == 3  % model set not complete: using three among five
    prompt = 'Model set not complete';
elseif scenario_index == 4  % only model uncertainty 
    prompt = 'Only uncertainty in model';
elseif scenario_index == 5  % model uncertainty + uncertain TPM
    prompt = 'Uncertainty in TPM and model';
elseif scenario_index == 6  % no uncertainty 
    prompt = 'No uncertainty ';
end
disp('+++++++++++++++++++++++++++++++++++++++++++++++++');
disp(['Scenario: ' num2str(scenario_index) ' -- ' prompt]);

%% Show Statistics
disp('+++++++++++++++++++++++++++++++++++++++++++++++++');
if length(err_dim) == 1
    if err_dim == 1
        disp('Error Statistics in Position');
    elseif err_dim == 2
        disp('Error Statistics in Speed');
    else
        error('main :: Unknown Error Dimensions (Check Point 1)');
    end
elseif length(err_dim) == 2
    disp('Error Statistics in Both Position and Speed');
else
    error('main :: Unknown Error Dimensions (Check Point 2)');
end

%% Display Error Plot
err_KF = zeros(1, EpisodeLength);
err_IMM_true = zeros(1, EpisodeLength);
err_IMM_nominal = zeros(1, EpisodeLength);
err_IMM_compensation = zeros(1, EpisodeLength);
err_IMM_Bayesian = zeros(1, EpisodeLength);
err_IMM_MaximumLikelihood = zeros(1, EpisodeLength);
err_IMM_RiskSensitive = zeros(1, EpisodeLength);
err_IMM_RobustKF = zeros(1, EpisodeLength);
err_DRIMM = zeros(1, EpisodeLength);
for episode = 1:NumMonteCarlo
    err_KF = err_KF + g_err_KF{episode};
    err_IMM_true = err_IMM_true + g_err_IMM_true{episode};
    err_IMM_nominal = err_IMM_nominal + g_err_IMM_nominal{episode};
    err_IMM_compensation = err_IMM_compensation + g_err_IMM_compensation{episode};
	err_IMM_Bayesian = err_IMM_Bayesian + g_err_IMM_Bayesian{episode};
    err_IMM_MaximumLikelihood = err_IMM_MaximumLikelihood + g_err_IMM_MaximumLikelihood{episode};
    err_IMM_RiskSensitive = err_IMM_RiskSensitive + g_err_IMM_RiskSensitive{episode};
    err_IMM_RobustKF = err_IMM_RobustKF + g_err_IMM_RobustKF{episode};
    err_DRIMM = err_DRIMM + g_err_DRIMM{episode};
end
err_KF = err_KF/NumMonteCarlo;
err_IMM_true = err_IMM_true/NumMonteCarlo;
err_IMM_nominal = err_IMM_nominal/NumMonteCarlo;
err_IMM_compensation = err_IMM_compensation/NumMonteCarlo;
err_IMM_Bayesian = err_IMM_Bayesian/NumMonteCarlo;
err_IMM_MaximumLikelihood = err_IMM_MaximumLikelihood/NumMonteCarlo;
err_IMM_RiskSensitive = err_IMM_RiskSensitive/NumMonteCarlo;
err_IMM_RobustKF = err_IMM_RobustKF/NumMonteCarlo;
err_DRIMM = err_DRIMM/NumMonteCarlo;

% Plot results
figure;
smt = 20; %50
% can also try "semilogx" instead of "plot"
plot(smooth(10*log10(err_KF), smt), 'k', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_IMM_true), smt), 'b', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_IMM_nominal), smt), 'g', 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_IMM_compensation), smt), 'm', 'LineWidth', 2); hold on;
if scenario_index ~= 4
plot(smooth(10*log10(err_IMM_Bayesian), smt), 'LineWidth', 2); hold on;
plot(smooth(10*log10(err_IMM_MaximumLikelihood), smt), 'LineWidth', 2); hold on;
end
plot(smooth(10*log10(err_IMM_RiskSensitive), smt), 'LineWidth', 2); hold on;
if scenario_index >= 3
plot(smooth(10*log10(err_IMM_RobustKF), smt), 'LineWidth', 2); hold on;
end
plot(smooth(10*log10(err_DRIMM), smt), 'r', 'LineWidth', 2); hold on;

% Plot Format
set(gca, 'FontSize', 18);
ylabel('Estimation Error (dB)','FontSize', 20, 'Interpreter', 'latex');
xlabel('Time Step','FontSize', 20, 'Interpreter', 'latex');
if scenario_index >= 3
    if scenario_index == 4
        leg1 = legend({'KF','IMM-T','IMM-N','IMM-C','IMM-RS','IMM-R','DRIMM'}, 'Interpreter', 'latex', 'Location', 'northwest');
    else
        leg1 = legend({'KF','IMM-T','IMM-N','IMM-C','IMM-B','IMM-M','IMM-RS','IMM-R','DRIMM'}, 'Interpreter', 'latex', 'Location', 'northwest');
    end
else
    leg1 = legend({'KF','IMM-T','IMM-N','IMM-C','IMM-B','IMM-M','IMM-RS','DRIMM'}, 'Interpreter', 'latex', 'Location', 'northwest');
end
set(leg1,'FontSize',15);

%% Display RMSE
RMSE_KF = 0;
RMSE_IMM_true = 0;
RMSE_IMM_nominal = 0;
RMSE_IMM_compensation = 0;
RMSE_IMM_Bayesian = 0;
RMSE_IMM_MaximumLikelihood = 0;
RMSE_IMM_RiskSensitive = 0;
RMSE_IMM_RobustKF = 0;
RMSE_DRIMM = 0;
for episode = 1:NumMonteCarlo
    RMSE_KF = RMSE_KF + g_RMSE_KF(episode);
    RMSE_IMM_true = RMSE_IMM_true + g_RMSE_IMM_true(episode);
    RMSE_IMM_nominal = RMSE_IMM_nominal + g_RMSE_IMM_nominal(episode);
    RMSE_IMM_compensation = RMSE_IMM_compensation + g_RMSE_IMM_compensation(episode);
	RMSE_IMM_Bayesian = RMSE_IMM_Bayesian + g_RMSE_IMM_Bayesian(episode);
    RMSE_IMM_MaximumLikelihood = RMSE_IMM_MaximumLikelihood + g_RMSE_IMM_MaximumLikelihood(episode);
    RMSE_IMM_RiskSensitive = RMSE_IMM_RiskSensitive + g_RMSE_IMM_RiskSensitive(episode);
    RMSE_IMM_RobustKF = RMSE_IMM_RobustKF + g_RMSE_IMM_RobustKF(episode);
    RMSE_DRIMM = RMSE_DRIMM + g_RMSE_DRIMM(episode);
end
RMSE_KF = RMSE_KF/NumMonteCarlo;
RMSE_IMM_true = RMSE_IMM_true/NumMonteCarlo;
RMSE_IMM_nominal = RMSE_IMM_nominal/NumMonteCarlo;
RMSE_IMM_compensation = RMSE_IMM_compensation/NumMonteCarlo;
RMSE_IMM_Bayesian = RMSE_IMM_Bayesian/NumMonteCarlo;
RMSE_IMM_MaximumLikelihood = RMSE_IMM_MaximumLikelihood/NumMonteCarlo;
RMSE_IMM_RiskSensitive = RMSE_IMM_RiskSensitive/NumMonteCarlo;
RMSE_IMM_RobustKF = RMSE_IMM_RobustKF/NumMonteCarlo;
RMSE_DRIMM = RMSE_DRIMM/NumMonteCarlo;
disp(['---------------------RMSE of KF: ' num2str(RMSE_KF)]);
disp(['---------------------RMSE of IMM-T (True): ' num2str(RMSE_IMM_true)]);
disp(['---------------------RMSE of IMM-N (Nominal): ' num2str(RMSE_IMM_nominal)]);
disp(['---------------------RMSE of IMM-C (Compensation): ' num2str(RMSE_IMM_compensation)]);
if scenario_index ~= 4
disp(['---------------------RMSE of IMM-B (Bayesian): ' num2str(RMSE_IMM_Bayesian)]);
disp(['---------------------RMSE of IMM-M (Maximum-Likelihood): ' num2str(RMSE_IMM_MaximumLikelihood)]);
end
disp(['---------------------RMSE of IMM-RS (Risk-sensitive): ' num2str(RMSE_IMM_RiskSensitive)]);
if scenario_index >= 3
disp(['---------------------RMSE of IMM-R (Robust): ' num2str(RMSE_IMM_RobustKF)]);
end
disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_DRIMM)]);


%% Display Time
time_KF = 0;
time_IMM_true = 0;
time_IMM_nominal = 0;
time_IMM_compensation = 0;
time_IMM_Bayesian = 0;
time_IMM_MaximumLikelihood = 0;
time_IMM_RiskSensitive = 0;
time_IMM_RobustKF = 0;
time_DRIMM = 0;
for episode = 1:NumMonteCarlo
    time_KF = time_KF + g_time_KF(episode);
    time_IMM_true = time_IMM_true + g_time_IMM_true(episode);
    time_IMM_nominal = time_IMM_nominal + g_time_IMM_nominal(episode);
    time_IMM_compensation = time_IMM_compensation + g_time_IMM_compensation(episode);
    time_IMM_Bayesian = time_IMM_Bayesian + g_time_IMM_Bayesian(episode);
    time_IMM_MaximumLikelihood = time_IMM_MaximumLikelihood + g_time_IMM_MaximumLikelihood(episode);
    time_IMM_RiskSensitive = time_IMM_RiskSensitive + g_time_IMM_RiskSensitive(episode);
    time_IMM_RobustKF = time_IMM_RobustKF + g_time_IMM_RobustKF(episode);
    time_DRIMM = time_DRIMM + g_time_DRIMM(episode);
end
time_KF = time_KF/NumMonteCarlo;
time_IMM_true = time_IMM_true/NumMonteCarlo;
time_IMM_nominal = time_IMM_nominal/NumMonteCarlo;
time_IMM_compensation = time_IMM_compensation/NumMonteCarlo;
time_IMM_Bayesian = time_IMM_Bayesian/NumMonteCarlo;
time_IMM_MaximumLikelihood = time_IMM_MaximumLikelihood/NumMonteCarlo;
time_IMM_RiskSensitive = time_IMM_RiskSensitive/NumMonteCarlo;
time_IMM_RobustKF = time_IMM_RobustKF/NumMonteCarlo;
time_DRIMM = time_DRIMM/NumMonteCarlo;
disp('+++++++++++++++++++++++++++++++++++++++++++++++++');
disp(['.....................Time of KF: ' num2str(time_KF)]);
disp(['.....................Time of IMM-T (True): ' num2str(time_IMM_true)]);
disp(['.....................Time of IMM-N (Nominal): ' num2str(time_IMM_nominal)]);
disp(['.....................Time of IMM-C (Compensation): ' num2str(time_IMM_compensation)]);
if scenario_index ~= 4
disp(['.....................Time of IMM-B (Bayesian): ' num2str(time_IMM_Bayesian)]);
disp(['.....................Time of IMM-M (Maximum-Likelihood): ' num2str(time_IMM_MaximumLikelihood)]);
end
disp(['.....................Time of IMM-RS (Risk-Sensitive): ' num2str(time_IMM_RiskSensitive)]);
if scenario_index >= 3
disp(['.....................Time of IMM-R (Robust): ' num2str(time_IMM_RobustKF)]);
end
disp(['.....................Time of DRIMM: ' num2str(time_DRIMM)]);


%% Generate Latex Tables
GenerateLatexTables;
