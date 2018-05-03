function C = loop_adapt_covariance(sigma, arx, arindex, xold, mu, c1, cmu, C, pc, hsig, cc, weights)
    artmp = (1/sigma) * (arx(:,arindex(1:mu))-repmat(xold,1,mu));
    
    C = (1-c1-cmu) * C ...                  % regard old matrix  
        + c1 * (pc*pc' ...                 % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
end