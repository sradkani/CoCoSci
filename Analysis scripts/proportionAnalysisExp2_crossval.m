function [coeffs,dev,stats] = proportionAnalysisExp2_crossval(alphabet, maxMotifLength, delta, alpha, change, remove1st, K)
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
    curves{i} = movmean(randomXCurve(...
        alphabet, maxMotifLength, sequences{i}, maxRandomX, delta, alpha), 3);
    
    diffcurves{i} = diff(curves{i});
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

for fold=1:K
    test_size = floor(length(sequences)/K);
    train_size = length(sequences) - test_size;
    train_inds = randperm(length(sequences), train_size);
    test_inds = 1:length(sequences);
    test_inds(sort(train_inds)) = [];
    
    train_sequences = sequences(train_inds);
    train_seqlengths = seqlengths(train_inds);
    test_sequences = sequences(test_inds);
    test_seqlengths = seqlengths(test_inds);
    train_curves = curves(train_inds);
    train_diffcurves = diffcurves(train_inds);
    test_curves = curves(test_inds);
    test_diffcurves = diffcurves(test_inds);
    
    % Make the variables for training
    % flatten diff curves and take absolute value
    train_diffs = abs(cell2mat(cellfun(@transpose, train_diffcurves, 'UniformOutput', false)));

    % get event positions
    train_eventpos = cell2mat(cellfun(@(x) (1:length(x))', train_diffcurves, 'UniformOutput', false));

    % get disengagement (1 = disengaged, 0 = continued)
    train_disengaged = ismember(1:sum(train_seqlengths), cumsum(train_seqlengths));

    % remove first two elements of disengaged (first gets lost when computing
    % curve, the second gets lost when computing the derivative)

    % disengagement except for last
    train_disengagementIdx = find(train_disengaged, sum(train_disengaged)-1);

    if change 
        train_removeIdx = [1:2, train_disengagementIdx+1, train_disengagementIdx+2];
    else
        train_removeIdx = [1, train_disengagementIdx+1];
    end

    train_disengaged(train_removeIdx) = [];

    % make disengaged column vector;
    train_disengaged = train_disengaged.';

    % remove 1st elements of diffs
    if remove1st
        idx = (eventpos==1 & diffs == 0);
        % remove diffs at eventpos 1
        train_diffs(idx) = [];
        train_disengaged(idx) = [];
        train_eventpos(idx) = [];

    end
    
    
    % Make the variables for testing
    % flatten diff curves and take absolute value
    test_diffs = abs(cell2mat(cellfun(@transpose, test_diffcurves, 'UniformOutput', false)));

    % get event positions
    test_eventpos = cell2mat(cellfun(@(x) (1:length(x))', test_diffcurves, 'UniformOutput', false));

    % get disengagement (1 = disengaged, 0 = continued)
    test_disengaged = ismember(1:sum(test_seqlengths), cumsum(test_seqlengths));

    % remove first two elements of disengaged (first gets lost when computing
    % curve, the second gets lost when computing the derivative)

    % disengagement except for last
    test_disengagementIdx = find(test_disengaged, sum(test_disengaged)-1);

    if change 
        test_removeIdx = [1:2, test_disengagementIdx+1, test_disengagementIdx+2];
    else
        test_removeIdx = [1, test_disengagementIdx+1];
    end

    test_disengaged(test_removeIdx) = [];

    % make disengaged column vector;
    test_disengaged = test_disengaged.';

    % remove 1st elements of diffs
    if remove1st
        idx = (eventpos==1 & diffs == 0);
        % remove diffs at eventpos 1
        test_diffs(idx) = [];
        test_disengaged(idx) = [];
        test_eventpos(idx) = [];

    end
    
    % fit the model to the training data
    [coeffs,dev,stats] = mnrfit([train_eventpos train_diffs], train_disengaged + 1);
    
    % find the AUC on the test data
    model_predictions = mnrval(coeffs, [test_eventpos test_diffs]);
    [~,~,~,auc] = perfcurve(test_disengaged + 1, model_predictions(:,2), 2);    
    AUCs(fold) = auc;
    
end
    
figure;
plot(AUCs);

mean_AUC = mean(AUCs)
