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

Data = importdata('Data.txt', '\t', 1);         % this is part of the data starting from 10s
                                                % use "Data (Speed-12).txt" instead if the drone's flying speed is 12m/s
EpisodeLength = length(Data.data(:,1));

%% Plot RTK Trajectory (Real Trajectory)
% plot3(ENU_GPS(1,:),ENU_GPS(2,:),ENU_GPS(3,:));
plot(Data.data(:,4), Data.data(:,5),'r','linewidth',2);
hold on;
plot(0,0,'b.','markersize',25,'MarkerFaceColor','b');
text(0,0,'Start');
plot(Data.data(end,4), Data.data(end,5),'s','markersize',8,'MarkerFaceColor','k');
text((Data.data(end,4)), (Data.data(end,5)),'End');
% legend('GPS','RTK');
set(gca, 'FontSize', 16);
ylabel('North (m)','FontSize', 18, 'Interpreter', 'latex');
xlabel('East (m)','FontSize', 18, 'Interpreter', 'latex');
axis equal;

%% Plot Real Velocity
figure;
% This is a data segment, starting from "10s"
plot(Data.data(:,1)+10, Data.data(:,6),'b','linewidth',2);
hold on;
plot(Data.data(:,1)+10, Data.data(:,7),'r-.','linewidth',2);
legend('East','North');
set(gca, 'FontSize', 16);
ylabel('Velocity (m/s)','FontSize', 18, 'Interpreter', 'latex');
xlabel('Time (s)','FontSize', 18, 'Interpreter', 'latex');
% figure;
% plot(Data.data(:,1), sqrt((Data.data(:,6)).^2 + (Data.data(:,7)).^2),'r-.','linewidth',2);

%% Load initialization
sys = LoadScenario_Real_Drone(EpisodeLength, Data.data);
x_hat0  = sys.init.x0;
P0      = sys.init.P0;
mu_hat0 = sys.init.mu0;

%% In which dimension do we calculate error statistics?
err_dim = 1:2;    % 1 (position), 2 (speed), 1:2 (all)

%% Standard IMM (Nominal Pi and Nominal a)
% IMM-N
tic
[x_hat_IMM_nominal, mu_hat_IMM_nominal] = IMM(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_nominal = toc/sys.K;
%disp(['Avg Time of IMM (Nominal): ' num2str(time_IMM_nominal)]);
err_Raw_IMM_nominal = sys.x(err_dim, :) - x_hat_IMM_nominal(err_dim, :);
err_IMM_nominal = sum(err_Raw_IMM_nominal.^2,1);
RMSE_IMM_nominal = sqrt(mean(err_IMM_nominal));
%disp(['---------------------RMSE of IMM (Nominal): ' num2str(RMSE_IMM_nominal)]);

%% Compensation IMM; Zhao Shunyi (2016) Journal of Franklin Institute
% IMM-C
% Recursive estimation for Markov jump linear systems with unknown transition probabilities: A compensation approach. 
tic
[x_hat_IMM_compensation, mu_hat_IMM_compensation] = IMM_Compensation(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_compensation = toc/sys.K;
%disp(['Avg Time of IMM (Compensation): ' num2str(time_IMM_compensation)]);
err_Raw_IMM_compensation = sys.x(err_dim, :) - x_hat_IMM_compensation(err_dim, :);
err_IMM_compensation = sum(err_Raw_IMM_compensation.^2,1);
RMSE_IMM_compensation = sqrt(mean(err_IMM_compensation));
%disp(['---------------------RMSE of IMM (Compensation): ' num2str(RMSE_IMM_compensation)]);

%% IMM with Bayesian estimation; Jilkov and Li (2004) IEEE Trans on Signal Processing
% IMM-B
% Online Bayesian estimation of transition probabilities for Markovian jump systems
tic
[x_hat_IMM_Bayesian, mu_hat_IMM_Bayesian] = IMM_Bayesian(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_Bayesian = toc/sys.K;
%disp(['Avg Time of IMM (Bayesian): ' num2str(time_IMM_Bayesian)]);
err_Raw_IMM_Bayesian = sys.x(err_dim, :) - x_hat_IMM_Bayesian(err_dim, :);
err_IMM_Bayesian = sum(err_Raw_IMM_Bayesian.^2,1);
RMSE_IMM_Bayesian = sqrt(mean(err_IMM_Bayesian));
%disp(['---------------------RMSE of IMM (Bayesian): ' num2str(RMSE_IMM_Bayesian)]);

%% IMM with Maximum Likelihood estimation; The Recursive-Kullback--Leibler (RKL) method; Orguner (2006)
% IMM-M
tic
[x_hat_IMM_MaximumLikelihood, mu_hat_IMM_MaximumLikelihood] = IMM_RKL(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_MaximumLikelihood = toc/sys.K;
%disp(['Avg Time of IMM (Maximum Likelihood, RKL): ' num2str(time_IMM_RKL)]);
err_Raw_MaximumLikelihood = sys.x(err_dim, :) - x_hat_IMM_MaximumLikelihood(err_dim, :);
err_IMM_MaximumLikelihood = sum(err_Raw_MaximumLikelihood.^2,1);
RMSE_IMM_MaximumLikelihood = sqrt(mean(err_IMM_MaximumLikelihood));
%disp(['---------------------RMSE of IMM (RKL): ' num2str(RMSE_IMM_RKL)]);

%% Risk-Sensitive IMM; Orguner (2008)
% IMM-RS
tic
[x_hat_IMM_RiskSensitive, mu_hat_IMM_RiskSensitive] = IMM_RS(sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_RiskSensitive = toc/sys.K;
%disp(['Avg Time of IMM (Risk-Sensitive): ' num2str(time_IMM_RS)]);
err_Raw_IMM_RiskSensitive = sys.x(err_dim, :) - x_hat_IMM_RiskSensitive(err_dim, :);
err_IMM_RiskSensitive = sum(err_Raw_IMM_RiskSensitive.^2,1);
RMSE_IMM_RiskSensitive = sqrt(mean(err_IMM_RiskSensitive));
%disp(['---------------------RMSE of IMM (RS): ' num2str(RMSE_IMM_RS)]);

%% IMM using Robust Kalman Filter; Wang (2021)
% IMM-R
tic
[x_hat_IMM_RobustKF, mu_hat_IMM_RobustKF] = DRIMM('P', sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_IMM_RobustKF = toc/sys.K;
%disp(['Avg Time of IMM (R): ' num2str(time_IMM_R)]);
err_Raw_IMM_RobustKF = sys.x(err_dim, :) - x_hat_IMM_RobustKF(err_dim, :);
err_IMM_RobustKF = sum(err_Raw_IMM_RobustKF.^2,1);
RMSE_IMM_RobustKF = sqrt(mean(err_IMM_RobustKF));
%disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_IMM_R)]);

%% Distributionally Robust IMM; This paper
% DRIMM
tic
% Possible Values of "Algo": 'mu', 'P', 'mu+P'
Algo = 'mu+P';
[x_hat_DRIMM, mu_hat_DRIMM] = DRIMM(Algo, sys.F, sys.G, sys.H, sys.Q, sys.R, sys.Pi, sys.y, sys.a, x_hat0, P0, mu_hat0);
time_DRIMM = toc/sys.K;
%disp(['Avg Time of DRIMM: ' num2str(time_DRIMM)]);
err_Raw_DRIMM = sys.x(err_dim, :) - x_hat_DRIMM(err_dim, :);
err_DRIMM = sum(err_Raw_DRIMM.^2,1);
RMSE_DRIMM = sqrt(mean(err_DRIMM));
%disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_DRIMM)]);

disp('+++++++++++++++++++++++++++++++++++++++++++++++++');
disp('Scenario: Real-World (Drone)');

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

%% Display RMSE
disp(['---------------------RMSE of IMM-N (Nominal): ' num2str(RMSE_IMM_nominal)]);
disp(['---------------------RMSE of IMM-C (Compensation): ' num2str(RMSE_IMM_compensation)]);
disp(['---------------------RMSE of IMM-B (Bayesian): ' num2str(RMSE_IMM_Bayesian)]);
disp(['---------------------RMSE of IMM-M (Maximum-Likelihood): ' num2str(RMSE_IMM_MaximumLikelihood)]);
disp(['---------------------RMSE of IMM-RS (Risk-sensitive): ' num2str(RMSE_IMM_RiskSensitive)]);
disp(['---------------------RMSE of IMM-R (Robust): ' num2str(RMSE_IMM_RobustKF)]);
disp(['---------------------RMSE of DRIMM: ' num2str(RMSE_DRIMM)]);


%% Display Time
disp('+++++++++++++++++++++++++++++++++++++++++++++++++');
disp(['.....................Time of IMM-N (Nominal): ' num2str(time_IMM_nominal)]);
disp(['.....................Time of IMM-C (Compensation): ' num2str(time_IMM_compensation)]);
disp(['.....................Time of IMM-B (Bayesian): ' num2str(time_IMM_Bayesian)]);
disp(['.....................Time of IMM-M (Maximum-Likelihood): ' num2str(time_IMM_MaximumLikelihood)]);
disp(['.....................Time of IMM-RS (Risk-Sensitive): ' num2str(time_IMM_RiskSensitive)]);
disp(['.....................Time of IMM-R (Robust): ' num2str(time_IMM_RobustKF)]);
disp(['.....................Time of DRIMM: ' num2str(time_DRIMM)]);


%% Generate Latex Tables
GenerateLatexTables;
