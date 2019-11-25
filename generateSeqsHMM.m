function [sequences, states, deltas, alphas] =...
    generateSeqsHMM(sizeAlphabet, maxMotifLength, seqLength, numSeqs)

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
        transitionMat(:,1+numPrevStates:max(pointers(:))) = C(seq)*alphas(seq).^i; 

        % insert 0's for mid motif states (they can only be reached through
        % within motif transitions)
        if i > 1 % because i=1 are single state motifs
            % get column indices of states that aren't beginnings of motifs
            nonBeginnerIdx = numPrevStates + elCountMat(:,2:end);
            transitionMat(:, nonBeginnerIdx(:)) = 0;
        end

        % convert to linear index to insert targeted deltas
        transitionMat(sub2ind(size(transitionMat),idx(:,1), idx(:,2))) = deltas(seq);

    end
    
    
    %% generate sequences
       
    % delta to C*alpha
    prior = [C(seq).*alphas(seq), transitionMat(1,2:end)];      
    
     % randomly sample one state between 1 and numStates with probabilities specified by prior
    startAt = randsample(numStates, 1, true, prior);

    % generate sequences
    [sequences(seq,:), states(seq,:)] = hmmgenerate2(seqLength, transitionMat, emissions, startAt,...
        'Symbols', cell2mat(cellfun(@(x) reshape(x.', 1, []), motifs, 'UniformOutput', false)));
    
end

csvwrite('states.csv', states)
csvwrite('sequences.csv', sequences)
csvwrite('deltas.csv', deltas)
csvwrite('alphas.csv', alphas)

end