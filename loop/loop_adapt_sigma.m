function sigma = loop_adapt_sigma(sigma, cs, damps, ps, chiN)
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); 
end