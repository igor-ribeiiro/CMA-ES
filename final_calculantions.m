function [eigeneval, B, C, D, invsqrtC, xmin, valmin] = final_calculantions(counteval, eigeneval, lambda, c1, cmu, N, C, D, numberOfLoops, arx, arindex, ...
                                                                        arfitness)
    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)'; % enforce symmetry
        [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
        D = sqrt(diag(D));        % D is a vector of standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
        fprintf("On iteration %d: \n", numberOfLoops);
        xmin = arx(:, arindex(1))' % Just printing the best result on each iteration
        valmin = arfitness(1)
    end
end