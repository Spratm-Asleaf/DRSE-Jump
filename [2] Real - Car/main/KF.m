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

function [x_hat, P] = KF(F,G,H,Q,R,Y,J,a,x_hat0,P0)
[~, episode_length] = size(Y);
P = P0;
X = x_hat0;
n = length(x_hat0);
x_hat = zeros(n, episode_length);

% filter starts
for k = 1:episode_length
    X = F*X + G*a(J(k));
    Z_ = H*X;

    P = F*P*F' + G*Q*G';

    S = H*P*H' + R;
    K = P*H'*(H*P*H' + R)^-1;
    X = X + K*(Y(k) - Z_);
    
    x_hat(:,k) = X;

    %P = (eye(n,n) - K*H)*P*(eye(n,n) - K*H)' + K*R*K';
    P = P - P*H'*S^(-1)*H*P;
end
