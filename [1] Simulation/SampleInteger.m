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

function index = SampleInteger(mu, count)
% Sample categorical indices according to the discrete distribution "mu".
% Return a column vector "index".
    N = length(mu);
    if abs(sum(mu)-1) > 1e-8 || max(mu) > 1 || min(mu) < 0
        error('Error in model probability vector!');
    end
    
    index = zeros(count, 1);
    
    for cnt = 1:count
        value = rand;
        summation = 0;
        for i = 1:N
            summation = summation + mu(i);
            if value <= summation
                index(cnt) = i;
                break;
            end
        end
    end
end
