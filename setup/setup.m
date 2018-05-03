function [lambda, mu, weights, mueff, cc, cs, c1, cmu, damps, pc, ps, B, D, C, invsqrtC, eigeneval, chiN, arx, arfitness] = setup(N)
  % Strategy parameter setting: Selection
  [lambda, mu, weights, mueff] = setup_selection(N);

  % Strategy parameter setting: Adaptation
  [cc, cs, c1, cmu, damps] = setup_adaptation(N, mueff);
  
  % Initialize dynamic (internal) strategy parameters and constants
  [pc, ps, B, D, C, invsqrtC, eigeneval, chiN, arx, arfitness] = setup_initialize_strategy_parameters(N, lambda);

end