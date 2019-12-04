function randomXs = findRandomXFitParams(alphabet, maxMotifLength, sequences, delta, alpha)
    
    maxRandomX = zeros(1,20);
    for seqLen=4:2:20
        maxRandomX(seqLen) = findMaxRandomX(alphabet, seqLen, maxMotifLength, delta, alpha);
    end
    
    randomXs = nan(1,length(sequences));
    % compute random(X) at each point in the sequence
    for i = 1:length(sequences)
        randomXs(i) = findRandomness(alphabet, maxMotifLength, sequences{i}, delta, alpha);

        % normalize by maxRandomX
        randomXs(i) = randomXs(i) ./ maxRandomX(length(sequences{i}));

    end
end