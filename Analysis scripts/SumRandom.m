function [coeffs,dev,stats, AUCs] = SumRandom(alphabet, maxMotifLength, delta, alpha, smooth, K)
    % this function takes the data for experiment 2 and visualizes the
    % relationship between proportion of disengaging and the ...
    datatable = parsejsPsychCSVExp2;
    sequences = datatable.sequences;

    seqlengths = cellfun(@length, sequences)';

    % remove sequences of length 1,2,29,30 
    sequences = sequences(seqlengths > 2 & seqlengths < 29, :);
    seqlengths = cellfun(@length, sequences).';

    if any(seqlengths == 30)
        error('there are sequences that were not terminated, write code to handle this before proceeding')    
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
        
        sumcurves{i} = cumsum(curves{i});
        
    end
    
    % cross validation loop
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
        train_sumcurves = sumcurves(train_inds);
        test_curves = curves(test_inds);
        test_sumcurves = sumcurves(test_inds);

        %% Make the variables for training
        % flatten diff curves and take absolute value
        train_SC = cell2mat(train_sumcurves)';

        % get event positions
        train_eventpos = cell2mat(cellfun(@(x) (1:length(x)), train_sumcurves, 'UniformOutput', false))';

        % get disengagement (1 = disengaged, 0 = continued)
        train_disengaged = ismember(1:sum(train_seqlengths), cumsum(train_seqlengths));

        % remove first two elements of disengaged (first gets lost when computing
        % curve, the second gets lost when computing the derivative)

        % disengagement except for last
        train_disengagementIdx = find(train_disengaged, sum(train_disengaged)-1);
    
        % remove elements to match diffs length
        train_removeIdx = [1, train_disengagementIdx+1];
        train_disengaged(train_removeIdx) = [];

        % make disengaged column vector;
        train_disengaged = train_disengaged.';
        

        
        %% Make the variables for testing
        % flatten diff curves and take absolute value
        test_SC = cell2mat(test_sumcurves)';

        % get event positions
        test_eventpos = cell2mat(cellfun(@(x) (1:length(x)), test_sumcurves, 'UniformOutput', false))';

        % get disengagement (1 = disengaged, 0 = continued)
        test_disengaged = ismember(1:sum(test_seqlengths), cumsum(test_seqlengths));

        % remove first two elements of disengaged (first gets lost when computing
        % curve, the second gets lost when computing the derivative)

        % disengagement except for last
        test_disengagementIdx = find(test_disengaged, sum(test_disengaged)-1);

        test_removeIdx = [1, test_disengagementIdx+1];
        test_disengaged(test_removeIdx) = [];

        % make disengaged column vector;
        test_disengaged = test_disengaged.';

%         %% Plots
%         % pool diff curves in bins

%         pooledg = quantile(train_SC,0:0.2:1);
%         linedges = linspace(min(train_SC), max(train_SC), 6);
% 
%         [bins, edges] = discretize(train_SC, pooledg);
% 
%         proportionDisengaged = accumarray(bins, train_disengaged, [], @mean);
%         stdErrorDisengaged = accumarray(bins, train_disengaged, [], @(x) std(x) ./ sqrt( length(x)));
% 
%         binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';
% 
%         % get trial position as a function of random x
%         figure;
%         plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3)
%         xlabel('random(X) sum', 'FontSize', 20)
%         ylabel('Proportion disenganged', 'FontSize', 20) 

        %% fit the model to the training data
        [coeffs,dev,stats] = mnrfit([train_eventpos train_SC train_SC.^2], train_disengaged + 1);

        %% find the AUC on the test data
        model_predictions = mnrval(coeffs, [test_eventpos test_SC test_SC.^2]);
        [~,~,~,auc] = perfcurve(test_disengaged + 1, model_predictions(:,2), 2);    
        AUCs(fold) = auc;
    end

    

end