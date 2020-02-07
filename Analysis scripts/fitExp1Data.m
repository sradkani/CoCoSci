%% Load and format human data
clear;clc;
datatable = parseJsPsychCSVExp1(1);

datatable(strcmp(datatable.stimulus, ''), :) = [];
datatable(strcmp(datatable.response, 'nan'), :) = [];
uniqueIDs = unique(datatable.trialID);
sequences = {};
subjRand_human = {};
for i=1:length(uniqueIDs)  
    seq = datatable.stimulus(datatable.trialID==uniqueIDs(i));
    sequences{i} = seq{1};
    subjRand_human{i} = datatable.response(datatable.trialID==uniqueIDs(i));
end

%% Set the model parameters
alphabet = 'ABC';
maxMotifLength = 3;

%% Find the optimal model parameters that fits the data the best
[modelParams, optCorr_each, optCorr_avg, modelHumanCorr_each, modelHumanCorr_avg] = ...
                        optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength);

%% Plot the heatmap
figure(1)
clims = [0.1, 0.7];
imagesc([0.25, 0.95], [0.05, 0.95], modelHumanCorr_each, clims)
colorbar
xlabel('delta');
ylabel('alpha');
title('model - human correlation')

figure(2)
clims = [0.4, 0.9];
imagesc([0.25, 0.95], [0.05, 0.95], modelHumanCorr_avg, clims)
h = colorbar;
ylabel(h, 'Correlation', 'FontSize', 26)
xlabel('delta', 'FontSize', 26);
ylabel('alpha', 'FontSize', 26);
ax = gca;
ax.FontSize = 20;
title('model - average human correlation', 'FontSize', 34);

save('plots/modelHumanCorr_each.mat', 'modelHumanCorr_each')
save('plots/modelHumanCorr_avg.mat', 'modelHumanCorr_avg')
