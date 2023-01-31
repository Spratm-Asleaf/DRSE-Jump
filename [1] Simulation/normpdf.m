function likelihood = normpdf(x, mu, cov)
    n = length(x);
    if isrow(x)
        x = x';
    end
    if isrow(mu)
        mu = mu';
    end
    likelihood = (2*pi)^(-n/2) * (det(cov))^(-1/2) * exp(-0.5 * (x - mu)'*cov^-1*(x - mu));
end