function xmin=cmaes   % (mu/mu_w, lambda)-CMA-ES
% --------------------  Initialization --------------------------------  
% User defined input parameters (need to be edited)
strfitnessfct = 'fitness_function';  % name of objective/fitness function
N = 10;               % number of objective variables/problem dimension
xmean = rand(N,1);    % objective variables initial point
sigma = 0.3;          % coordinate wise standard deviation (step size)
stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
stopeval = 1e3*N^2;   % stop after stopeval number of function evaluations

% Initial setup for the algorithm
addpath(genpath('setup'));
addpath(genpath('loop'));
[lambda, mu, weights, mueff, cc, cs, c1, cmu, damps, pc, ps, B, D, C, invsqrtC, eigeneval, chiN, arx, arfitness] = setup(N);

% -------------------- Generation Loop --------------------------------
% the next 40 lines contain the 20 lines of interesting code 
counteval = 0; % Is the number of function evaluations. Is limited by stopeval
numberOfLoops = 0; % How many loops on the script

while counteval < stopeval
    numberOfLoops = numberOfLoops + 1;
    
    % All the loop steps that I don't understand
    [arx, arfitness, counteval, arindex, xmean, ps, pc, C, sigma] = loop_proccessing(lambda, xmean, sigma, B, C, D, N, strfitnessfct, ...
                                                                               counteval, arx, arfitness, weights, mu, cs, ps, mueff, invsqrtC, ...
                                                                               chiN, cc, pc, ...
                                                                               damps, c1, cmu);

    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    [eigeneval, B, C, D, invsqrtC, xmin, valmin] = final_calculantions(counteval, eigeneval, lambda, c1, cmu, N, C, D, numberOfLoops, arx, ...
                                                                         arindex, arfitness);

    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
    if valmin <= stopfitness || max(D) > 1e7 * min(D)
        break;
    end

end % while, end generation loop

xmin = xmean'; % Return best point of last iteration.
               % Notice that xmean is expected to be even
               % better. It was xmin before, I did put the xmean because he
               % says it is better.
                             
end %function
