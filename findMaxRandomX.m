function [maxRandomX, maxSeq] = findMaxRandomX(alphabet, seqLength, maxMotifLength, delta, alpha)
    % this function finds the sequence built with a given set of alphabet 
    % with length seqLength that has the maximum subjective
    % randomness under a specific parsing model with parameters
    % maxMotifLength, delta and alpha
        
    sizeAlphabet = length(alphabet);
    
    % for small alphabet size and sequence length we can generate all
    % possible sequences, otherwise we will sample from all possible
    % sequences that are generated from a fair alphabetSize coin!
    if sizeAlphabet^seqLength < 100000
        sequences = permn(alphabet, seqLength);
    else
        idx = randi(sizeAlphabet, 50000, seqLength);
        sequences = arrayfun(@(x) alphabet(x), idx);
    end
    
    random_Xs = findRandomness(alphabet, maxMotifLength, sequences, delta, alpha);
    [maxRandomX, maxIdx] = max(random_Xs);
    maxSeq = sequences(maxIdx,:);
    
end