function [arx, arfitness, counteval, arindex, xmean, ps, pc, C, sigma] = loop_proccessing(lambda, xmean, sigma, B, C, D, N, strfitnessfct, counteval, arx, arfitness, weights, mu, cs, ps, mueff, invsqrtC, chiN, cc, pc, damps, c1, cmu)
    % Generate and evaluate lambda offspring
    [arx, arfitness, counteval] = loop_evaluate_function(lambda, xmean, sigma, B, D, N, strfitnessfct, counteval, arx, arfitness);

    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex, xold, xmean] = loop_sort_and_weight(arfitness, xmean, arx, weights, mu);

    % Cumulation: Update evolution paths
    [ps, hsig, pc] = loop_update_evolution_paths(N, cs, ps, mueff, invsqrtC, xmean, xold, sigma, counteval, lambda, chiN, cc, pc);

    % Adapt covariance matrix C
    C = loop_adapt_covariance(sigma, arx, arindex, xold, mu, c1, cmu, C, pc, hsig, cc, weights);

    % Adapt step size sigma
    sigma = loop_adapt_sigma(sigma, cs, damps, ps, chiN);
end