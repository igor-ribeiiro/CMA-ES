function [ps, hsig, pc] = loop_update_evolution_paths(N, cs, ps, mueff, invsqrtC, xmean, xold, sigma, counteval, lambda, chiN, cc, pc)
      ps = (1-cs)*ps ... 
            + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma; 
      hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
      pc = (1-cc)*pc ...
            + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
end