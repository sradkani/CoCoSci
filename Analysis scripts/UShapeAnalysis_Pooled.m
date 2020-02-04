%% Define the parameters
    clear;clc;
    alphabet = ['ABC'];
    maxMotifLength = 3;
    delta = 0.75;
    alpha = 0.35;

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

%% find model randomness estimates
    % get randomness curves
    curves = cell(size(sequences, 1), 1);

    for i = 1:size(sequences,1)
        curves{i} = randomXCurve(...
            alphabet, maxMotifLength, sequences{i}, maxRandomX, delta, alpha);
    end

    % flatten randomness curves
    randoms = cell2mat(cellfun(@transpose, curves, 'UniformOutput', false));
    
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

%% Bin the randomness values
    % pool randomness values in bins
    pooledg = quantile(randoms,[0 0.2 0.4 0.6 0.8 1]);
    linedges = linspace(min(randoms), max(randoms), 6);

    [bins, edges] = discretize(randoms, pooledg);

    proportionDisengaged = accumarray(bins, disengaged, [], @mean);
    stdErrorDisengaged = accumarray(bins, disengaged, [], @(x) std(x) ./ sqrt( length(x)));

    binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

    % plot proportionDisengaged as a function of model randomness estimates
    figure;
    plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3)
    xlabel('Random(x)', 'FontSize', 20)
    ylabel('Proportion disenganged', 'FontSize', 20)



























