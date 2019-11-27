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
    % compute denominator for C (sum over number of motifs times alpha^length)
    denom = arrayfun(@(x, y) x .* alphas.^y,...
        motifSizes(:,2).', 1:length(motifs), 'UniformOutput', false);

    % (1-deltas)/ 
    % (3 .* alpha.^1 + 6 * alpha.^2 + ... + numOfMaxLengthMotifs * alpha^maxMotifLength)
    C = (1-deltas) ./ sum(cell2mat(denom'), 1);

    
    transitionMat = generateTransitionMat(...
        motifs, motifSizes, numStates, maxMotifLength, alpha, delta, C);

    prior = [C.*alpha, transitionMat(1,2:end)];     
    TRANS_HAT = [0 prior; zeros(size(transitionMat,1),1) transitionMat];
    EMIS_HAT = [zeros(1,size(emissions,2)); emissions];
    
    %% compute random(x)
   for seq=1:size(sequences,1) 
        States = hmmviterbi(sequences(seq,:), TRANS_HAT, EMIS_HAT, ...
                              'Symbols', alphabet);                  
        States = States - 1;
        % find P(X|regular)
        idx = buffer(States, 2, 1, 'nodelay');
        P_X_regular = prior(States(1));
        P_X_regular = P_X_regular * prod( transitionMat(sub2ind(size(transitionMat), idx(1,:), idx(2,:))), 2);
        
        % the base of logarithm scales the randomness measure
        % for now use natural logarithm
        l_X = log(length(sequences(seq,:)));
        random_X(seq) = -l_X - log(P_X_regular);
   end
end













