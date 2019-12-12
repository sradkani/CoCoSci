function [modelParams, optCorr_each, optCorr_avg, modelHumanCorr_each, modelHumanCorr_avg] = optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength, delta0, alpha0)
    
    modelParams.alphabet = alphabet;
    modelParams.maxMotifLength = maxMotifLength;
    
%     fun = @(params) -1 * findModelHumanCorr(sequences, subjRand_human, modelParams, params(1), params(2));
%     [optParams, optCorr] = fminsearch(fun, [delta0, alpha0]);
    
%     optCorr = -optCorr;
%     modelParams.delta = optParams(1);
%     modelParams.alpha = optParams(2);

    alpha = 0.05:0.1:0.95;
    delta = 0.25:0.1:0.95;
    for i=1:length(alpha)
        for j=1:length(delta)
            [corr_each, corr_avg] = findModelHumanCorr(sequences, subjRand_human, modelParams, delta(j), alpha(i));
            modelHumanCorr_each(i,j) = corr_each;
            modelHumanCorr_avg(i,j) = corr_avg;
        end
    end
    
    [optCorr_each, maxInd_each] = max(modelHumanCorr_each(:));
    [I,J] = ind2sub(size(modelHumanCorr_each), maxInd_each);
    modelParams.delta_each = delta(J);
    modelParams.alpha_each = alpha(I);
    
    [optCorr_avg, maxInd_avg] = max(modelHumanCorr_avg(:));
    [I,J] = ind2sub(size(modelHumanCorr_avg), maxInd_avg);
    modelParams.delta_avg = delta(J);
    modelParams.alpha_avg = alpha(I);

end