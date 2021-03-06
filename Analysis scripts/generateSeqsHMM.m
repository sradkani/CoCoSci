function [sequences, states, deltas, alphas] =...
    generateSeqsHMM(alphabet, maxMotifLength, seqLength, numSeqs, paramRange, flatPrior)

%% generate motifs

motifs = generateMotifs(alphabet, maxMotifLength);

% get dimensions of each entry in motifs
motifSizes = cell2mat(cellfun(@(x) size(x), motifs, 'UniformOutput', false)');

% number of states (multiply dimensions along columns and sum)
numStates = sum(prod(motifSizes, 2));


%% set parameters for sequence generation sequences with different parameter settings
% (delta = staying in motif, alpha = prior over motifs);

% give extremal dispersion of parameters
if paramRange
    
    if  numSeqs <= 1
        error('can only give range of parameters if numSeqs > 1')
  
    elseif mod(numSeqs,2) ~= 0 
        error('numSeqs must be even if more than one sequenc length is specified')
       
    end
    
    % two gaussians 
    deltas = [0.2+0.6*randn(1, numSeqs/2), (0.2*randn(1, numSeqs/2) + 0.95)];
    
    % set values below 0 and above 1 to 0.01 and 0.99
    deltas(deltas <= 0.01) = 0.01;
    deltas(deltas >= 0.99) = 0.99;
    
    alphas = rand(1,numSeqs);
else
    deltas = rand(1,numSeqs); alphas = rand(1, numSeqs);
end


% compute denominators for C (sum over number of motifs times alpha^length)
motifs = generateMotifs('ABC', 3);

% get number of non-zero entries for each motif size
nonZero = cellfun(@(x) size(x, 1), motifs);

% duplicate and deduct identity (because whenever we are in a motif of
% length i, delta will replace one of the C*a^i)
coeffs = repmat(nonZero, [maxMotifLength 1]) - eye(maxMotifLength);

% multiply with alpha^i where i = motiflength
denom = coeffs * alphas.^((1:maxMotifLength).');

% generate C (will have dimensions maxMotifLength x 1)
C = (1-deltas) ./ denom;

% repeat (first 3 rows will have first entry of C, 12 will have second
% etc.)
C = repelem(C, cellfun(@(x) numel(x), motifs), 1);

% emissions are deterministic
emissions = eye(numStates);

% preallocate sequence and state storage
if length(seqLength) > 1
    sequences = cell(numSeqs*length(seqLength), 1);
    states = cell(numSeqs*length(seqLength), 1);
else
    sequences = nan(numSeqs, seqLength);
    states = nan(numSeqs, seqLength);

end


for seq = 1:numSeqs
    %% generate transition matrix
 
    transitionMat = generateTransitionMat(...
    motifs, motifSizes, numStates, maxMotifLength, alphas(seq), deltas(seq), C(:, seq));
    
    %% generate sequences
        
    if flatPrior
        % flat prior assigns equal probability acorss motifs (normalized by
        % number of motifs)
        prior = ones(1, size(transitionMat(1,:), 1)) ./...
            repelem(prod(motifSizes,2), prod(motifSizes,2));
    else
        % delta to C*alpha
        prior = [C(seq).*alphas(seq), transitionMat(1,2:end)];
    end
    
    % loop through several sequence lengths
    if length(seqLength) > 1

        for i = 1:length(seqLength)
            
         % randomly sample one state between 1 and numStates with probabilities specified by prior
         startAt = randsample(numStates, 1, true, prior);
       
        [sequences{(seq-1)*length(seqLength) + i}, states{(seq-1)*length(seqLength) + i}] =...
            hmmgenerate2(seqLength(i), transitionMat, emissions, startAt,...
            'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
        
         % redo if there is a duplicate
        while (seq-1)*length(seqLength) + i >...
                size(uniqueRowsCA(sequences(1:(seq-1)*length(seqLength) + i), 'rows'),1)
            
        [sequences{(seq-1)*length(seqLength) + i}, states{(seq-1)*length(seqLength) + i}] =...
            hmmgenerate2(seqLength(i), transitionMat, emissions, startAt,...
            'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
        
         % disp('replaced a duplicate')
        end
        end
        
    else
    
         % randomly sample one state between 1 and numStates with probabilities specified by prior
        startAt = randsample(numStates, 1, true, prior);

        [sequences(seq,:), states(seq,:)] =...
            hmmgenerate2(seqLength, transitionMat, emissions, startAt,...
            'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
        
        % redo if there is a duplicate
        while seq > size(unique(sequences(1:seq, :), 'rows'),1)
             [sequences(seq,:), states(seq,:)] =...
            hmmgenerate2(seqLength, transitionMat, emissions, startAt,...
            'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
        
         % disp('replaced a duplicate')
        end
        
    end
end

% remove duplicate values
if ischar(alphabet)
    sequences = char(sequences);
end

% shuffle rows
sequences = sequences(randperm(size(sequences, 1)),: );

%csvwrite('states.csv', states)
csvwrite('/Users/galraz1/Developer/CoCoSci/Experiment2/sequencesExp2.csv', sequences)
%csvwrite('deltas.csv', deltas)
%csvwrite('alphas.csv', alphas)

end