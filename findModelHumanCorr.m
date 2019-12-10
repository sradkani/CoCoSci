function ModelHumanCorrelation = findModelHumanCorr(sequences, subjRand_human, modelParams, delta, alpha)
    % This function finds the correlation between model randomness estimates and
    % human subjective randomness for a given set of model parameters
    % each sequence has a couple of subjRAnd_human associated with it
    % subjRand_human is a cell array, each cell containing human data
    % correpsonding to each sequence
    
    maxRandomX = zeros(1,20);
    for seqLen=4:2:20
        maxRandomX(seqLen) = findMaxRandomX(modelParams.alphabet, seqLen, modelParams.maxMotifLength, delta, alpha);
    end
    
    subjRand_model = nan(length(sequences), 1);
    % compute random(X) at each point in the sequence
    for i = 1:length(sequences)
        subjRand_model(i,1) = findRandomness(modelParams.alphabet, modelParams.maxMotifLength, sequences{i}, delta, alpha);

        % normalize by maxRandomX
        subjRand_model(i,1) = subjRand_model(i,1) ./ maxRandomX(length(sequences{i}));
        
    end
    
    % unroll the human data and repeat subjRand_model accordingly
%     subjRand_human_unrolled = cell2mat(cellfun(@(x) reshape(x.', 1, []), subjRand_human, 'UniformOutput', false))';
%     subjRand_model_repeated = [];
%     for i=1:length(subjRand_model)
%         subjRand_model_repeated = [subjRand_model_repeated; repmat(subjRand_model(i), [length(subjRand_human{i}), 1])];
%     end
    
    % find the correlation between human and model subjective randomness 
%     CorrMat = corr([subjRand_human_unrolled, subjRand_model_repeated]);
    
    % take the average of human ratings for each sequence
    subjRand_human_mean = cellfun(@(x) mean(x), subjRand_human);
    
    CorrMat = corr([subjRand_human_mean', subjRand_model]);
    ModelHumanCorrelation = CorrMat(1,2);
    
    display(alpha)
    display(delta)
    display(ModelHumanCorrelation)
end



