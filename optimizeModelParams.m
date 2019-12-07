function [modelParams, optCorr] = optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength, delta0, alpha0)
    
    modelParams.alphabet = alphabet;
    modelParams.maxMotifLength = maxMotifLength;
    
    fun = @(params) -1 * findModelHumanCorr(sequences, subjRand_human, modelParams, params(1), params(2));
    [optParams, optCorr] = fminsearch(fun, [delta0, alpha0]);
    
    optCorr = -optCorr;
    modelParams.delta = optParams(1);
    modelParams.alpha = optParams(2);
end