function xmin=cmaes   % (mu/mu_w, lambda)-CMA-ES
  % --------------------  Initialization --------------------------------  
  % User defined input parameters (need to be edited)
  strfitnessfct = 'frosenbrock';  % name of objective/fitness function
  N = 10;               % number of objective variables/problem dimension
  xmean = rand(N,1);    % objective variables initial point
  sigma = 0.3;          % coordinate wise standard deviation (step size)
  stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
  stopeval = 1e3*N^2;   % stop after stopeval number of function evaluations
  
  
  % Initial setup for the algorithm
  addpath(genpath('setup'));
  [lambda, mu, weights, mueff, cc, cs, c1, cmu, damps, pc, ps, B, D, C, invsqrtC, eigeneval, chiN] = setup(N);

  % -------------------- Generation Loop --------------------------------
  % the next 40 lines contain the 20 lines of interesting code 
  counteval = 0; % Is the number of function evaluations. Is limited by stopeval
  numberOfLoops = 0; % How many loops on the script
  addpath(genpath('loop'))
  while counteval < stopeval
      numberOfLoops = numberOfLoops + 1;
      % Generate and evaluate lambda offspring
      for k=1:lambda,
          arx(:,k) = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C) 
          arfitness(k) = feval(strfitnessfct, arx(:,k)); % objective function call
          counteval = counteval+1;
      end
      
      % Sort by fitness and compute weighted mean into xmean
      [arfitness, arindex, xold, xmean] = loop_sort_and_weight(arfitness, xmean, arx, weights, mu);
      
    
      % Cumulation: Update evolution paths
      [ps, hsig, pc] = loop_update_evolution_paths(N, cs, ps, mueff, invsqrtC, xmean, xold, sigma, counteval, lambda, chiN, cc, pc);
      

      % Adapt covariance matrix C
      artmp = (1/sigma) * (arx(:,arindex(1:mu))-repmat(xold,1,mu));
      C = (1-c1-cmu) * C ...                  % regard old matrix  
           + c1 * (pc*pc' ...                 % plus rank one update
                   + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
           + cmu * artmp * diag(weights) * artmp'; % plus rank mu update

      % Adapt step size sigma
      sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
    
      % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
      if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
          eigeneval = counteval;
          C = triu(C) + triu(C,1)'; % enforce symmetry
          [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
          D = sqrt(diag(D));        % D is a vector of standard deviations now
          invsqrtC = B * diag(D.^-1) * B';
          fprintf("On iteration %d: \n", numberOfLoops);
          xmin = arx(:, arindex(1))' % Just printing the best result on each iteration
      end
    
      % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
      if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
          break;
      end

  end % while, end generation loop

  xmin = arx(:, arindex(1))'; % Return best point of last iteration.
                             % Notice that xmean is expected to be even
                             % better.
% ---------------------------------------------------------------  
function f=frosenbrock(x)
    if size(x,1) < 2 error('dimension must be greater one'); end
    f = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);