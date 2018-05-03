function [arx, arfitness, counteval] = loop_evaluate_function(lambda, xmean, sigma, B, D, N, strfitnessfct, counteval, arx, arfitness)
    for k=1:lambda
        arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
        arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call
        counteval = counteval+1;
    end
end