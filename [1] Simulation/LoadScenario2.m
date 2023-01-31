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

function [sys] = LoadScenario2(EpisodeLength)
    K = EpisodeLength;  % time steps of each episode
    sys.K = K;

    n = 2;      % dimension of state vector x
    m = 1;      % dimension of measurement vector y
    T = 1;      % sampling time
    N = 4;      % number of candidate models
    sys.n = n; sys.m = m; sys.T = T; sys.N = N;

    x = zeros(n, K);    % state vector
    y = zeros(m, K);    % measurement vector
    j = zeros(1, K);    % model index

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
    sys.Q = 2^2;
    sys.R = 100^2;

    % system init state
    x0 = [80000; 400];
    init.x0 = x0;

    % system init state covariance
    init.P0 = [
        100^2 0
        0     100^2
    ];

    % system init model probability 
    mu0 = [0.6; 0.1; 0.1; 0.1;0.1];
    init.mu0 = mu0;

    % model transition probability matrix
    Pi = [
        0.6     0.1     0.1     0.1     0.1
        0.1     0.6     0.1     0.1     0.1
        0.1     0.1     0.6     0.1     0.1
        0.1     0.1     0.1     0.6     0.1
        0.1     0.1     0.1     0.1     0.6
    ];
    sys.Pi = Pi;

    % system first model
    j0 = SampleInteger(mu0, 1);
    init.j0 = j0;

    % save system init attibutes
    sys.init = init;

    % three acceleration modes
    a = [0;20;-20;40;-40];
    sys.a = a;

    % generate states and measurements
    for k = 1:K
        % system dynamics
        w = sqrt(sys.Q)*randn;
        if k == 1
            j(k) = SampleInteger(Pi(j0,:), 1);
            x(:, k) = sys.F * x0 + sys.G*(a(j(k)) + w);
        else
            j(k) = SampleInteger(Pi(j(k-1),:), 1);
            x(:, k) = sys.F * x(:, k-1) + sys.G*(a(j(k)) + w);
        end

        % measurement dynamics
        v = sqrt(sys.R)*randn;
        y(:, k) = sys.H * x(:, k) + v;
    end

    % save syste states and measurements
    sys.x = x;
    sys.y = y;
    sys.j = j;
end











