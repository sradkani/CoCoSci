function [randomnessMeasures] = findRandomness(sizeAlphabet, maxMotifLength, sequences, delta, alpha)

    %% generate motifs

    % initialize alphabet
    alphabet = 1:sizeAlphabet;

    % preallocate motif cell
    motifs = cell(1, maxMotifLength);

    % generate motifs
    for i = 1:maxMotifLength
        % each entry in motifs will be a matrix containing all motifs of a
        % given length
        motifs{i} = permn(alphabet, i);

        % remove motifs that can be described by shorter motifs

        % this is the case when longer motifs repeat themselves
        % e.g. ACAC can be represented by 'AC' motif alone

        % find partitions of motifs of length i 
        % (e.g. i = 4 -> mod(4, 0:3) -> 4 0 0 1 -> 
        % find 0's to get valid partitions)
        partitions = find(mod(i, 1:i-1) == 0); 

        % iterate through partition, e.g. if i = 4, partitions are 1 and 2
        for partition = partitions
            % duplicate 1:partition columns of matrix to check 
            % which rows contain duplications
            duplicatedMat = repmat(motifs{i}(:,1:partition), [1 i/partition]);

            % rows that are in duplicated matrix are redundant
            [~, redundantRows] = intersect(motifs{i}, duplicatedMat, 'rows');

            % delete redundant rows
            motifs{i}(redundantRows, :) = [];

        end
    end
    
    % get dimensions of each entry in motifs
    motifSizes = cell2mat(cellfun(@(x) size(x), motifs, 'UniformOutput', false)');

    % number of states (multiply dimensions along columns and sum)
    numStates = sum(prod(motifSizes, 2));
    
    % =============================== this should be modified
    C = (1-delta) / (2 * (alpha + alpha^2));        

    % emissions are deterministic
    emissions = zeros(numStates, sizeAlphabet);
    stateOutcomes = cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false));
    for i=1:sizeAlphabet
        emissions(:,i) = (stateOutcomes == alphabet(i))';
    end
    
    
    %% generate transition matrix
    
    % preallocate transition matrix
    transitionMat = nan(numStates);
    
    for i = 1:maxMotifLength

    % first put deltas in right place for each state (one per row)
    % the delta should in the place pointing to the state representing the next
    % state in the motif (or start the motif again)

        % get element count matrix (labels all elements 1 to numel(motifs{i}))
        elCountMat = reshape(1:numel(motifs{i}), flip(motifSizes(i,:))).';

        % 1 + to indicate that i move to the next state in the motif
        % flip and .' to reshape in row major order (instead of default col)
        pointers = 1 + elCountMat;

        % subtract the number of columns in motifSizes for last column to
        % re-point to the beginning of the motif
        pointers(:, end) = pointers(:,end) - motifSizes(i, 2);

        % add number of states from previous motif lengths
        numPrevStates = sum(prod(motifSizes(1:i-1,:), 2));
        pointers = pointers + numPrevStates;

          % get index for all within motif states
        idx = [(1+numPrevStates:max(pointers(:))).', reshape(pointers.', 1, []).'];

        % populate columns with C*a.^(length of motif this state is a part of), 
        % some will later be replaced with 0 or delta
        transitionMat(:,1+numPrevStates:max(pointers(:))) = C * alpha^i; 

        % insert 0's for mid motif states (they can only be reached through
        % within motif transitions)
        if i > 1 % because i=1 are single state motifs
            % get column indices of states that aren't beginnings of motifs
            nonBeginnerIdx = numPrevStates + elCountMat(:,2:end);
            transitionMat(:, nonBeginnerIdx(:)) = 0;
        end

        % convert to linear index to insert targeted deltas
        transitionMat(sub2ind(size(transitionMat),idx(:,1), idx(:,2))) = delta;

    end
    
    prior = [C.*alpha, transitionMat(1,2:end)];     
    TRANS_HAT = [0 prior; zeros(size(transitionMat,1),1) transitionMat];
    EMIS_HAT = [zeros(1,size(emissions,2)); emissions];
    
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
        randomnessMeasures(seq) = -l_X - log(P_X_regular);
   end
end













