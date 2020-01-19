function saveMaxRandomX(alphabet, maxMotifLength, maxSeqLength, delta, alpha)
% this creates a mat file with maxRandom(x) for a specified alpha and delta
% from length 2 to maxSeqLength -1

maxRandomX = nan(1, maxSeqLength-1);

% normalize by maxRandomX
for i = 2:30
    maxRandomX(i-1) = findMaxRandomX(alphabet, i, maxMotifLength, delta, alpha);
end

save(sprintf('savedMaxRandomX/maxRandomX_delta%.2f_alpha%.2f.mat', delta, alpha), 'maxRandomX')