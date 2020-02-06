function [coeffs,dev,stats,auc] = Exp2Analysis_randomX_unnormalized(alphabet, maxMotifLength, delta, alpha, remove1st, smooth)
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

for i = 1:size(sequences,1)
    if smooth
        curves{i} = movmean(randomXCurve_unnormalized(...
            alphabet, maxMotifLength, sequences{i}, delta, alpha), 3);
    else
        curves{i} = randomXCurve_unnormalized(...
            alphabet, maxMotifLength, sequences{i}, delta, alpha);
    end
    
end


% flatten diff curves and take absolute value
randomx = cell2mat(cellfun(@transpose, curves, 'UniformOutput', false));

% get event positions
eventpos = cell2mat(cellfun(@(x) (1:length(x))', curves, 'UniformOutput', false));

% get disengagement (1 = disengaged, 0 = continued)
disengaged = ismember(1:sum(seqlengths), cumsum(seqlengths));

% disengagement except for last
disengagementIdx = find(disengaged, sum(disengaged)-1);

% remove first element of disengaged gets lost when computing
removeIdx = [1, disengagementIdx+1];
disengaged(removeIdx) = [];

% make disengaged column vector;
disengaged = disengaged.';

% remove 1st elements of diffs
if remove1st
    idx = (eventpos==1 & randomx == 0);
    % remove diffs at eventpos 1
    randomx(idx) = [];
    disengaged(idx) = [];
    eventpos(idx) = [];
    
end

% pool diff curves in bins
pooledg = quantile(randomx,0:0.2:1);
linedges = linspace(min(randomx), max(randomx), 6);

[bins, edges] = discretize(randomx, pooledg);

proportionDisengaged = accumarray(bins, disengaged, [], @mean);
stdErrorDisengaged = accumarray(bins, disengaged, [], @(x) std(x) ./ sqrt( length(x)));

binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

% binning random(x) values and get mean proportion of disengaged
figure;
plot(randomx, eventpos,  'ro', 'LineWidth', 3)
xlabel('random(x)', 'FontSize', 20)
ylabel('length', 'FontSize', 20)


% get trial position as a function of random x
figure;
plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3)
xlabel('random(x)', 'FontSize', 20)
ylabel('Proportion disenganged', 'FontSize', 20) 

[coeffs,dev,stats] = mnrfit([eventpos randomx], disengaged + 1);
model_predictions = mnrval(coeffs, [eventpos randomx]);
[~,~,~,auc] = perfcurve(disengaged + 1, model_predictions(:,2), 2);  

plottable = table();

plottable.disengaged = disengaged;
plottable.diffs = randomx;

plottable.binMeans = binMeans(bins);

plottable.eventpos = eventpos;

writetable(plottable, 'Rdata2.csv')

