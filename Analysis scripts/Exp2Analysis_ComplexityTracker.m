function [coeffs,dev,stats, AUCs] = Exp2Analysis_ComplexityTracker(alphabet, maxMotifLength, delta, alpha, smooth)
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
    
        %% Make the variables for training
        % flatten diff curves and take absolute value
         CT = cell2mat(complexity_tracker)';

        % get event positions
         eventpos = cell2mat(cellfun(@(x) (1:length(x)),  complexity_tracker, 'UniformOutput', false))';

        % get disengagement (1 = disengaged, 0 = continued)
         disengaged = ismember(1:sum( seqlengths), cumsum( seqlengths));

        % remove first two elements of disengaged (first gets lost when computing
        % curve, the second gets lost when computing the derivative)

        % disengagement except for last
         disengagementIdx = find( disengaged, sum( disengaged)-1);
    
        % remove elements to match diffs length
         removeIdx = [1:2,  disengagementIdx+1,  disengagementIdx+2];
         disengaged( removeIdx) = [];

        % make disengaged column vector;
         disengaged =  disengaged.';
         
         % pool diff curves in bins
        pooledg = quantile(CT,0:0.2:1);
        linedges = linspace(min(CT), max(CT), 6);

        [bins, edges] = discretize(CT, pooledg);

        proportionDisengaged = accumarray(bins, disengaged, [], @mean);
        stdErrorDisengaged = accumarray(bins, disengaged, [], @(x) std(x) ./ sqrt( length(x)));

        binMeans = (0.5 * (edges(1:end-1) + edges(2:end)))';

        % binning random(x) values and get mean proportion of disengaged
        figure;
        plot(CT, eventpos,  'ro', 'LineWidth', 3)
        xlabel('CT value', 'FontSize', 20)
        ylabel('event position', 'FontSize', 20)


        % get trial position as a function of random x
        figure;
        plot(binMeans, proportionDisengaged,  'ro-', 'LineWidth', 3)
        xlabel('CT value', 'FontSize', 20)
        ylabel('Proportion disenganged', 'FontSize', 20) 


        %% fit the model to the training data
        [coeffs,dev,stats] = mnrfit([ eventpos  CT],  disengaged + 1);

        plottable = table();

        plottable.disengaged = disengaged;
        plottable.CT = CT;

        plottable.binMeans = binMeans(bins);

        plottable.eventpos = eventpos;

        writetable(plottable, 'Rdata_CT.csv')


    end

    

