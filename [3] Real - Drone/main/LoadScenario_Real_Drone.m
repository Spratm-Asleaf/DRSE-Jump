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

function [sys] = LoadScenario_Real_Drone(EpisodeLength, Data)
    K = EpisodeLength;  % time steps of each episode
    sys.K = K;
    
    n = 2;                          % dimension of state vector x
    m = 1;                          % dimension of measurement vector y
    T = 0.1;%Data(2,1) - Data(1,1);      % sampling time
    
    sys.n = n; sys.m = m; sys.T = T;

    % system matrices
    F = [
        1 T
        0 1
    ];
    sys.F = F;
    G = [
        T^2/2
        T
    ];
    sys.G = G;
    sys.H = [1 0];
    sys.Q = 0.01^2;
    sys.R = 0.1^2;

    % system init state
    x0 = [0,0]';
    init.x0 = x0;

    % system init state covariance
    init.P0 = [
        100^2   0 
        0       100^2 
    ];

    N = 3;                          % number of candidate models
    sys.N = N;
    % system init model probability 
    init.mu0 = [0.6; 0.2; 0.2];

    % model transition probability matrix
    sys.Pi = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i == j
                sys.Pi(i,j) = 0.8;
            else
                sys.Pi(i,j) = 0.1;
            end
        end
    end

    % save system init attibutes
    sys.init = init;

    % three acceleration modes
    sys.a = [0   2.5  -2.5];

    % load syste states and measurements
    rng(46);
    aixs = 'e';
    if aixs == 'e'          % east
        sys.x = (Data(:, [4,6]))';
        sys.y = (Data(:, 2))';              % there is no built-in pre-denoising in this GPS solution
    elseif aixs == 'n'      % north
        sys.x = (Data(:, [5,7]))';
        sys.y = (Data(:, 3))';              % there is no built-in pre-denoising in this GPS solution
    else
        error('LoadScenario_Real_Drone :: Error in axis!');
    end
end











