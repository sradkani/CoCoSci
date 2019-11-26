function  transitionMat = generateTransitionMat(...
    motifs, motifSizes, numStates, maxMotifLength, alpha, delta, C)

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
    