function [coeffs,dev,stats,auc] = proportionAnalysisExp2(alphabet, maxMotifLength, delta, alpha, change, remove1st, smooth)
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
        curves{i} = movmean(randomXCurve(...
            alphabet, maxMotifLength, sequences{i}, delta, alpha), 3);
    else
        curves{i} = randomXCurve(...
            alphabet, maxMotifLength, sequences{i}, delta, alpha);
    end
    
    diffcurves{i} = diff(curves{i});
    
    % find the complexity tracker
    slope_tracker = diffcurves{i} > 0;
    slope_change = [0, diff(slope_tracker)];
    idx = [find(slope_change), length(slope_change)];
    vec = 1:idx(1);
    for j=2:length(idx)
        vec = [vec, 1:(idx(j)-idx(j-1))];
    end
    complexity_tracker{i} = vec;
    end

% % remove the seqs with change in randomness < 0.04 in eventpos<4
% k = 1;
% for i=1:length(sequences)
%     if length(diffcurves{i}) >1
%         if any(abs(diffcurves{i}(1:2))<0.04)
%             continue;
%         else
%             temp{k,1} = sequences{i};
%             temp{k,2} = diffcurves{i};
%             temp{k,3} = seqlengths(i);
%             k = k + 1;
%         end
%     end
% end
% sequences = temp(:,1);
% diffcurves = temp(:,2);
% seqlengths = cell2mat(temp(:,3));

if change
    0;
else
    diffcurves = curves;
end

% flatten diff curves and take absolute value
diffs = cell2mat(cellfun(@transpose, diffcurves, 'UniformOutput', false));
        
% complexity tracker
CT = abs(cell2mat(complexity_tracker))';

% get event positions
eventpos = cell2mat(cellfun(@(x) (1:length(x))', diffcurves, 'UniformOutput', false));

% get disengagement (1 = disengaged, 0 = continued)
disengaged = ismember(1:sum(seqlengths), cumsum(seqlengths));

% remove first two elements of disengaged (first gets lost when computing
% curve, the second gets lost when computing the derivative)

% disengagement except for last
disengagementIdx = find(disengaged, sum(disengaged)-1);

if change 
    removeIdx = [1:2, disengagementIdx+1, disengagementIdx+2];
else
    removeIdx = [1, disengagementIdx+1];
end

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
pooledg = quantile(diffs,0:0.2:1);
linedges = linspace(min(diffs), max(diffs), 6);

[bins, edges] = discretize(diffs, pooledg);

proportionDisengaged = accumarray(bins, disengaged, [], @mean);
stdErrorDisengaged = accumarray(bins, disengaged, [], @(x) std(x) ./ sqrt( length(x)));

binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

% binning random(x) values and get mean proportion of disengaged
figure;
plot(diffs, eventpos,  'ro', 'LineWidth', 3)
if change
    xlabel('Change in random(x)', 'FontSize', 20)
else
    xlabel('random(x)', 'FontSize', 20)
end

ylabel('length', 'FontSize', 20)


% get trial position as a function of random x
figure;
plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3)
if change
    xlabel('complexity tracker value', 'FontSize', 20)
else
    xlabel('random(x)', 'FontSize', 20)
end

ylabel('Proportion disenganged', 'FontSize', 20) 

[coeffs,dev,stats] = mnrfit([eventpos diffs], disengaged + 1);
model_predictions = mnrval(coeffs, [eventpos diffs]);
[~,~,~,auc] = perfcurve(disengaged + 1, model_predictions(:,2), 2);  

plottable = table();

plottable.disengaged = disengaged;
plottable.diffs = diffs;

plottable.binMeans = binMeans(bins);

plottable.eventpos = eventpos;

writetable(plottable, 'Rdata2.csv')
% 
% plottable = table();
% 
% plottable.binMeans = binMeans;
% plottable.propDisengaged = proportionDisengaged;
% plottable.stdErrorDisengaged = stdErrorDisengaged;
% 
% writetable(plottable, 'Rdata.csv')

% 
