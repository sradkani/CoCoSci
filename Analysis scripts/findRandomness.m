function [random_X] = findRandomness(alphabet, maxMotifLength, sequences, delta, alpha)

% this function computes random(seq) for sequences in rows of 'sequences'

    %% check data integrity
    
    sizeAlphabet = length(alphabet);
    
    % number of unique values in each sequence
    numUniq = cellfun(@(c) length(unique(c)), num2cell(sequences, 2));
    
    % if any of the sequences have more unique values than the size of the
    % alphabet, throw error
    if any(numUniq > sizeAlphabet)
        error(strcat('Sequence(s) #', sprintf(' %d', find(numUniq > sizeAlphabet)),...
        ' have more unique values than the specified size of the alphabet'))
    end

    %% generate motifs

    motifs = generateMotifs(alphabet, maxMotifLength);
    
    % get dimensions of each entry in motifs
    motifSizes = cell2mat(cellfun(@(x) size(x), motifs, 'UniformOutput', false)');

    % number of states (multiply dimensions along columns and sum)
    numStates = sum(prod(motifSizes, 2));
    

    %% generate emissions
    
    emissions = zeros(numStates, sizeAlphabet);
    stateOutcomes = cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false));
    for i=1:sizeAlphabet
        emissions(:,i) = (stateOutcomes == alphabet(i))';
    end
    
%% generate transition matrix

    % compute denominators for C (sum over number of motifs times alpha^length)
    motifs = generateMotifs('ABC', 3);
    
    % get number of non-zero entries for each motif size
    nonZero = cellfun(@(x) size(x, 1), motifs);
    
    % duplicate and deduct identity (because whenever we are in a motif of
    % length i, delta will replace one of the C*a^i)
    coeffs = repmat(nonZero, [maxMotifLength 1]) - eye(maxMotifLength);
    
    % multiply with alpha^i where i = motiflength
    denom = coeffs * alpha.^(1:maxMotifLength).';
    
    % generate C (will have dimensions maxMotifLength x 1)
    C = (1-delta) ./ denom;
    
    % repeat (first 3 rows will have first entry of C, 12 will have second
    % etc.)
    C = repelem(C, cellfun(@(x) numel(x), motifs));
    
    transitionMat = generateTransitionMat(...
        motifs, motifSizes, numStates, maxMotifLength, alpha, delta, C);

    prior = [C(1).*alpha, transitionMat(1,2:end)];     
    TRANS_HAT = [0 prior; zeros(size(transitionMat,1),1) transitionMat];
    EMIS_HAT = [zeros(1,size(emissions,2)); emissions];
    
    %% compute random(x)
   for seq=1:size(sequences,1) 
        States = hmmviterbi(sequences(seq,:), TRANS_HAT, EMIS_HAT, ...
                              'Symbols', alphabet);                  
        States = States - 1;
        % find P(X|regular)
        idx = buffer(States, 2, 1, 'nodelay');
        P_X_regular = prior(States(1)) *...
            prod( transitionMat(sub2ind(size(transitionMat), idx(1,:), idx(2,:))), 2);
        
        % the base of logarithm scales the randomness measure
        % for now use natural logarithm
        l_X = length(sequences(seq,:));
        random_X(seq) = -l_X - log(P_X_regular)/log(3);
   end
end