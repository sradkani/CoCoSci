% this function generates sequences from an HMM with different transition
% probabilities
function seqMat = generateSeqs(sizeAlphabet, motifLength, seqLength, numSeqs)

% initialize alphabet
alphabet = 1:sizeAlphabet;

% preallocate motif cell
motifs = cell(1, motifLength);

% motifs of length motifLength to 1
for i = 1:motifLength
    % each entry in motifs will be a matrix containing all motifs of a
    % given length
    motifs{i} = permn(alphabet, i);
end

% # of motifs for each motif length
motifNumPerLength = sizeAlphabet.^(1:motifLength);

% total number of motifs
% numIndices = sum(motifNumPerLength);

% preallocate matrix for sequences
seqMat = nan(numSeqs, seqLength);

% range of transition probabilities from 0 to 1
transitionProbs = linspace(0, 1, numSeqs);

for i = 1:numSeqs
    
    seq = [];
    
    % sample transition probability for this sequence
    TP = transitionProbs(i);

    % initialize at random motif
    % linIdx = randi([1 numIndices]);
    
    % initialize at random motif length and random motif
    motifLen = randi([1 sizeAlphabet]);
    motifNum = randi([1 length(motifs{motifLen})]);
    
    while length(seq) <= seqLength
        
        % with sample new index with probability TP
        if rand < TP
            
            motifLen = randi([1 sizeAlphabet]);
            motifNum = randi([1 length(motifs{motifLen})]);
        % get linear index to for transition
        %    linIdx = randi([1 numIndices]);
        end
        
        % check which motif length to sample from
        % motifLen = sum(linIdx > cumsum(motifNumPerLength)) + 1;

        % check which motif to sample from from this motif length
        % motifNum = linIdx - sum(motifNumPerLength(1:motifLen-1));

        % index into motifs
        currentMotif = motifs{motifLen}(motifNum, :);

        % concatenate current motif with sequence
        seq = [seq currentMotif];
    end
    
    % shorten sequence if its longer than sequence length
    seq = seq(1:seqLength);
    
    seqMat(i,:) = seq;
end


csvwrite('sequences.csv', seqMat)
