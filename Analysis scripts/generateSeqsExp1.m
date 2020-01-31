alphabet = 'ABC';
maxMotifLength = 3;
numSeqs = 14;
paramRange = 1;
flatPrior = 1;

seqlengths = 4:2:20;

sequences = string();

for i = 1:length(seqlengths)
    sequences = strvcat(sequences,...
        generateSeqsHMM(alphabet, maxMotifLength, seqlengths(i), numSeqs, paramRange, flatPrior));
end

csvwrite('/Users/galraz1/Developer/CoCoSci/Experiment1/sequencesExp1.csv', sequences)
