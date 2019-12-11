function curves = getCurvesForExp2Seqs(alphabet, maxMotifLength, delta, alpha)
% gets curves for exp2 sequences
% tries to find file in folder, if it's not there, it computes it itself

sequences = cell2mat(cellfun(@(x) strrep(x,',',''),...
    importdata('sequencesExp2.csv'), 'UniformOutput', false));

try
maxRandomX = importdata(sprintf('maxRandomX_delta%.2f_alpha%.2f.mat', delta, alpha));
catch
    warning('maxRandomX not stored for these parameter settings. Computing them now, which might take longer')
    maxRandomX = nan(1, size(sequences, 2));

    % normalize by maxRandomX
    for i = 2:size(sequences, 2)
        maxRandomX(i-1) = findMaxRandomX(alphabet, i, maxMotifLength, delta, alpha);
    end
end

curves = nan(size(sequences, 1), size(sequences, 2)-1);

for i = 1:size(sequences,1)
    curves(i,:) = randomXCurve(...
        alphabet, maxMotifLength, sequences(i,:), maxRandomX, delta, alpha);

end


curveTable = table();
curveTable.curve1 = curves(9, :)';
curveTable.sequence1 = sequences(9, :)';

curveTable.curve2 = curves(14, :)';
curveTable.sequence2 = sequences(9, :)';

curveTable.curve3 = curves(34, :)';
curveTable.sequence3 = sequences(9, :)';

curveTable.t = (1:size(curves, 2))';

writetable(curveTable, 'curveRData.csv');


figure;
for i = 1:25
    subplot(5,5,i)
    plot(1:size(curves, 2), curves(i,:), 'LineWidth', 2)
    title(sequences(i,:), 'FontSize', 9)
    ylim([0 1])
end


figure;
for i = 1:25
    subplot(5,5,i)
    plot(1:size(curves, 2)-1, diff(curves(i,:)), 'LineWidth', 2)
    title(sequences(i,:), 'FontSize', 9)
end