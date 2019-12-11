function [modelParams, optCorr, modelHumanCorr] = optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength, delta0, alpha0)
    
    modelParams.alphabet = alphabet;
    modelParams.maxMotifLength = maxMotifLength;
    
%     fun = @(params) -1 * findModelHumanCorr(sequences, subjRand_human, modelParams, params(1), params(2));
%     [optParams, optCorr] = fminsearch(fun, [delta0, alpha0]);
    
%     optCorr = -optCorr;
%     modelParams.delta = optParams(1);
%     modelParams.alpha = optParams(2);

    alpha = 0.05:0.05:0.95;
    delta = 0.25:0.05:0.95;
    for i=1:length(alpha)
        for j=1:length(delta)
            modelHumanCorr(i,j) = findModelHumanCorr(sequences, subjRand_human, modelParams, delta(j), alpha(i));
        end
    end
    
    [optCorr, maxInd] = max(modelHumanCorr(:));
    [I,J] = ind2sub(size(modelHumanCorr), maxInd);
    modelParams.delta = delta(J);
    modelParams.alpha = alpha(I);

end