function curves = getCurvesForExp2Seqs(alphabet, maxMotifLength, sequences, delta, alpha)

maxRandomX = nan(1, size(sequences, 2));

% normalize by maxRandomX
for i = 2:size(sequences, 2)
    maxRandomX(i-1) = findMaxRandomX(alphabet, i, maxMotifLength, delta, alpha);
end

curves = nan(size(sequences, 1), size(sequences, 2)-1);

for i = 1:size(sequences,1)
    curves(i,:) = randomXCurve(...
        alphabet, maxMotifLength, sequences(i,:), maxRandomX, delta, alpha);

end

figure;
for i = 1:25
    subplot(5,5,i)
    title(sequences(i,:), 'FontSize', 6)
    plot(1:size(curves, 2), curves(i,:), 'LineWidth', 2)
end