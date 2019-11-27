%% Load the sequences
sequences = csvread('sequences.csv');

%% Define the parameters
sizeAlphabet = 3;
alphabet = 1:sizeAlphabet;
maxMotifLength = 4;
delta = 0.5;
alpha = 0.1;

%% Find randomness curve for each sequence
numSeq = size(sequences, 1);  seqLength = size(sequences, 2);
randomXs = nan(numSeq, seqLength);

for i=1:numSeq
    randomXs(i,:) = randomXCurve(alphabet, maxMotifLength, sequences(i), delta, alpha);
end

%% Save randomness curves
csvwrite('randomXs.csv', randomXs);