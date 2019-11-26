function [sequences, states, deltas, alphas] =...
    generateSeqsHMM(alphabet, maxMotifLength, seqLength, numSeqs)

%% generate motifs

motifs = generateMotifs(alphabet, maxMotifLength);

% get dimensions of each entry in motifs
motifSizes = cell2mat(cellfun(@(x) size(x), motifs, 'UniformOutput', false)');

% number of states (multiply dimensions along columns and sum)
numStates = sum(prod(motifSizes, 2));


%% set parameters for sequence generation sequences with different parameter settings
% (delta = staying in motif, alpha = prior over motifs);
deltas = rand(1,numSeqs); alphas = rand(1, numSeqs);
C = (1-deltas) ./ (2.* (alphas + alphas.^2));

% emissions are deterministic
emissions = eye(numStates);

% preallocate sequence and state storage
sequences = nan(numSeqs, seqLength);
states = nan(numSeqs, seqLength);

for seq = 1:numSeqs
    %% generate transition matrix
 
    transitionMat = generateTransitionMat(...
    motifs, motifSizes, numStates, maxMotifLength, alphas(seq), deltas(seq), C(seq));
    
    %% generate sequences
       
    % delta to C*alpha
    prior = [C(seq).*alphas(seq), transitionMat(1,2:end)];
    
     % randomly sample one state between 1 and numStates with probabilities specified by prior
    startAt = randsample(numStates, 1, true, prior);

    % generate sequences
    [sequences(seq,:), states(seq,:)] =...
        hmmgenerate2(seqLength, transitionMat, emissions, startAt,...
        'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
    
end

csvwrite('states.csv', states)
csvwrite('sequences.csv', sequences)
csvwrite('deltas.csv', deltas)
csvwrite('alphas.csv', alphas)

end