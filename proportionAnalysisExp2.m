function proportionAnalysisExp2(alphabet, maxMotifLength, delta, alpha, nBins)
% this function takes the data for experiment 2 and visualizes the
% relationship between proportion of disengaging and the diff of random x
datatable = parsejsPsychCSVExp2;
sequences = datatable.sequences;
seqlengths = cellfun(@length, sequences)';

if any(seqlengths == 30)
    error('there are sequences that were not terminated, write code to handle this before proceeding')
end

% tries to find file containing maxRandomX with these parameters in folder, 
% if it's not there, it computes it itself
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

% get curves and derivatives
curves = cell(size(sequences, 1), 1);
diffcurves =  cell(size(sequences, 1), 1);

for i = 1:size(sequences,1)
    curves{i} = randomXCurve(...
        alphabet, maxMotifLength, sequences{i}, maxRandomX, delta, alpha);
    
    diffcurves{i} = diff(curves{i});
end

% flatten diff curves and take absolute value
diffs = abs(cell2mat(cellfun(@transpose, diffcurves, 'UniformOutput', false)));

% get disengagement (1 = disengaged, 0 = continued)
disengaged = ismember(1:sum(seqlengths), cumsum(seqlengths));

% remove first two elements of disengaged (first gets lost when computing
% curve, the second gets lost when computing the derivative)

% disengagement except for last
disengagementIdx = find(disengaged, sum(disengaged) - 1);

removeIdx = [1:2, disengagementIdx+1, disengagementIdx+2];

disengaged(removeIdx) = [];

% pool diff curves in bins
[bins, edges] = discretize(diffs, linspace(min(diffs), max(diffs), nBins));

proportionDisengaged = accumarray(bins, disengaged', [], @mean);

binMeans = 0.5 * (edges(1:end-1) + edges(2:end));

figure;
plot(binMeans, proportionDisengaged, 'r*')
xlabel('Change in random(x)')
ylabel('Proportion disenganged')



