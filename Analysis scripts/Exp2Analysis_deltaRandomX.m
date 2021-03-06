function [coeffs,dev,stats,auc] = Exp2Analysis_deltaRandomX(alphabet, maxMotifLength, delta, alpha, remove1st, smooth)
% this function takes the data for experiment 2 and visualizes the
% relationship between proportion of disengaging and the diff of random x
datatable = parsejsPsychCSVExp2;
sequences = datatable.sequences;

seqlengths = cellfun(@length, sequences)';

% remove sequences of length 1,2,29,30 
sequences = sequences(seqlengths > 2 & seqlengths < 29, :);
seqlengths = cellfun(@length, sequences).';

if any(seqlengths == 30)
    error('there are sequences that were not terminated, write code to handle this before proceeding')    
end

% tries to find file containing maxRandomX with these parameters in folder, 
% if it's not there, it computes it itself
try
maxRandomX = importdata(sprintf('savedMaxRandomX/maxRandomX_delta%.2f_alpha%.2f.mat', delta, alpha));
catch
    warning('maxRandomX not stored for these parameter settings. Computing them now, which might take longer')
    maxRandomX = nan(1, size(sequences, 2));

    % normalize by maxRandomX
    for i = 2:max(cellfun(@length, sequences))
        maxRandomX(i-1) = findMaxRandomX(alphabet, i, maxMotifLength, delta, alpha);
    end
end

% get curves and derivatives
curves = cell(size(sequences, 1), 1);
diffcurves =  cell(size(sequences, 1), 1);

for i = 1:size(sequences,1)
    if smooth
        curves{i} = movmean(randomXCurve_normalized(...
            alphabet, maxMotifLength, sequences{i}, maxRandomX, delta, alpha), 3);
    else
        curves{i} = randomXCurve_normalized(...
            alphabet, maxMotifLength, sequences{i}, maxRandomX, delta, alpha);
    end
    
    diffcurves{i} = diff(curves{i});

end
  

% flatten diff curves and take absolute value
diffs = abs(cell2mat(cellfun(@transpose, diffcurves, 'UniformOutput', false)));
        
% get event positions
eventpos = cell2mat(cellfun(@(x) (1:length(x))', diffcurves, 'UniformOutput', false));

% get disengagement (1 = disengaged, 0 = continued)
disengaged = ismember(1:sum(seqlengths), cumsum(seqlengths));

% remove first two elements of disengaged (first gets lost when computing
% curve, the second gets lost when computing the derivative)

% disengagement except for last
disengagementIdx = find(disengaged, sum(disengaged)-1);

removeIdx = [1:2, disengagementIdx+1, disengagementIdx+2];

disengaged(removeIdx) = [];

% make disengaged column vector;
disengaged = disengaged.';

% remove 1st elements of diffs
if remove1st
    idx = (eventpos==1 & diffs == 0);
    % remove diffs at eventpos 1
    diffs(idx) = [];
    disengaged(idx) = [];
    eventpos(idx) = [];
    
end

% pool diff curves in bins
pooledg = quantile(diffs,0:0.1:1);
linedges = linspace(min(diffs), max(diffs), 6);

[bins, edges] = discretize(diffs, pooledg);

proportionDisengaged = accumarray(bins, disengaged, [], @mean);
stdErrorDisengaged = accumarray(bins, disengaged, [], @(x) std(x) ./ sqrt(length(x)));

binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

% binning random(x) values and get mean proportion of disengaged
figure;
plot(diffs, eventpos,  'ro', 'LineWidth', 4)
set(gca, 'FontSize', 12); box off
xlabel('Change in random(x)', 'FontSize', 20)
ylabel('event position', 'FontSize', 20)


% get trial position as a function of random x
figure;
errorbar(binMeans, proportionDisengaged, stdErrorDisengaged, 'bo-', 'LineWidth', 3)
set(gca, 'FontSize', 12); box off
xlabel('Change in random(x)', 'FontSize', 20)
ylabel('Proportion disenganged', 'FontSize', 20) 

[coeffs,dev,stats] = mnrfit([eventpos diffs], disengaged + 1);
model_predictions = mnrval(coeffs, [eventpos diffs]);
[~,~,~,auc] = perfcurve(disengaged + 1, model_predictions(:,2), 2);  

plottable = table();

plottable.disengaged = disengaged;
plottable.diffs = diffs;

plottable.binMeans = binMeans(bins);
plottable.err = stdErrorDisengaged(bins);

plottable.eventpos = eventpos;

writetable(plottable, 'Rdata_deltaRandomX.csv')
