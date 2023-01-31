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

function [sys] = LoadScenario_Real_Car(EpisodeLength, Data)
    K = EpisodeLength;  % time steps of each episode
    sys.K = K;
    
    n = 2;                          % dimension of state vector x
    m = 1;                          % dimension of measurement vector y
    T = 0.1;%Data(2,1) - Data(1,1);      % sampling time
    
    sys.n = n; sys.m = m; sys.T = T;

    % system matrices
    sys.F = [ 
        1 T
        0 1
    ];
    sys.G = [
        T^2/2
        T
    ];
    sys.H = [1 0];
    sys.Q = 0.1^2;
    sys.R = 1^2;

    % system init state
    x0 = (Data(1, 2:3))';
    init.x0 = x0;

    % system init state covariance
    init.P0 = [
        100^2 0
        0     100^2
    ];

    N = 3;                          % number of candidate models
    sys.N = N;
    % system init model probability 
    mu0 = [0.8; 0.1; 0.1];
    init.mu0 = mu0;

    % model transition probability matrix
    sys.Pi = [
        0.8     0.1     0.1
        0.1     0.8     0.1
        0.1     0.1     0.8
    ];

    % save system init attibutes
    sys.init = init;

    % three acceleration modes
    a = [0;10;-10];
    sys.a = a;

    % load syste states and measurements
    rng(46);
    aixs = 'e';
    if aixs == 'e'          % east
        sys.x = (Data(:, [4,6]))';
        sys.y = (Data(:, 2) + randn(EpisodeLength,1))';
        % N.B.: Since mature commercial GPS solutions have build-in de-noising functions,
        %       some noises are added back here to recover very original GPS positions.
    elseif aixs == 'n'      % north
        sys.x = (Data(:, [5,7]))';
        sys.y = (Data(:, 3) + randn(EpisodeLength,1))';
    else
        error('LoadScenario_Real_1 :: Error in axis!');
    end
end











