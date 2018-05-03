function [arfitness, arindex, xold, xmean] = loop_sort_and_weight(arfitness, xmean, arx, weights, mu)
      [arfitness, arindex] = sort(arfitness); % minimization
      xold = xmean;
      xmean = arx(:,arindex(1:mu))*weights;   % recombination, new mean value
end