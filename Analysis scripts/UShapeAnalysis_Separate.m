%% Define the parameters
    clear;clc;
    alphabet = ['ABC'];
    maxMotifLength = 3;
    delta = 0.2;
    alpha = 0.5;

%% load the data
    datatable = parsejsPsychCSVExp2;
    sequences = datatable.sequences;
    seqlengths = cellfun(@length, sequences)';

    % remove sequences of length 1 
    sequences = sequences(seqlengths > 2, :);
    seqlengths = cellfun(@length, sequences).';

    if any(seqlengths == 30)
        error('there are sequences that were not terminated, write code to handle this before proceeding')
    end

    % tries to find file containing maxRandomX with these parameters in folder, 
    % if it's not there, it computes it itself
%     try
%     maxRandomX = importdata(sprintf('savedMaxRandomX/maxRandomX_delta%.2f_alpha%.2f.mat', delta, alpha));
%     catch
%         warning('maxRandomX not stored for these parameter settings. Computing them now, which might take longer')
%         maxRandomX = nan(1, size(sequences, 2));
% 
%         % normalize by maxRandomX
%         for i = 2:max(cellfun(@length, sequences))
%             maxRandomX(i-1) = findMaxRandomX(alphabet, i, maxMotifLength, delta, alpha);
%         end
%     end

%% find model randomness estimates
    % get randomness curves
    curves = cell(size(sequences, 1), 1);

    for i = 1:size(sequences,1)
        curves{i} = randomXCurve(...
            alphabet, maxMotifLength, sequences{i}, delta, alpha);
    end

    % flatten randomness curves
    randoms = cell2mat(cellfun(@transpose, curves, 'UniformOutput', false));
    
    % make the eventpositions
    eventpositions = cell2mat(arrayfun(@(x) [2:x], seqlengths, 'UniformOutput', false));
    
%% find human disengagement behavior
    % get disengagement (1 = disengaged, 0 = continued)
    disengaged = ismember(1:sum(seqlengths), cumsum(seqlengths));

    % remove first element of disengaged (gets lost when computing curve)
    % disengagement except for last
    disengagementIdx = find(disengaged, sum(disengaged) - 1);
    removeIdx = [1, disengagementIdx+1];
    disengaged(removeIdx) = [];

    % make disengaged column vector;
    disengaged = disengaged.';

%% Group randomness values and disengagement behavior based on event position
    % group the randomness values for same eventpositions
    randomness_values = {};
    disengagements = {};
    for i=1:max(eventpositions)
        randomness_values{i} = randoms(eventpositions==i);
        disengagements{i} = disengaged(eventpositions==i);
    end
   
    
%% Plot proportionDisengaged as a funciton of random(X) for each seqlength separately
    figure(2);
    % evenposition=1 is empty
    for i=2:length(randomness_values)
        % pool randomness values in bins
        pooledg = quantile(randomness_values{i},[0 0.2 0.4 0.6 0.8 1]);
        linedges = linspace(min(randomness_values{i}), max(randomness_values{i}), 6);

        [bins, edges] = discretize(randomness_values{i}, pooledg);

        proportionDisengaged = accumarray(bins, disengagements{i}, [], @mean);
        stdErrorDisengaged = accumarray(bins, disengagements{i}, [], @(x) std(x) ./ sqrt( length(x)));

        binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

        % plot proportionDisengaged as a function of model randomness estimates
        subplot(6,5,i-1);
        plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3);
        title(strcat('seqLength', num2str(i)))
%         xlabel('Random(x)', 'FontSize', 20)
%         ylabel('Proportion disenganged', 'FontSize', 20)
    end
    



























