%% Load and format human data
clear;clc;
parseJsPsychCSVExp1

datatable(strcmp(datatable.stimulus, ''), :) = [];
datatable(strcmp(datatable.response, 'nan'), :) = [];
uniqueIDs = unique(datatable.trialID);
sequences = {};
subjRand_human = {};
for i=1:length(uniqueIDs)  
    seq = datatable.stimulus(datatable.trialID==uniqueIDs(i));
    sequences{i} = seq{1};
    subjRand_human{i} = str2num(cell2mat(datatable.response(datatable.trialID==uniqueIDs(i))));
end



%% Set the model parameters
alphabet = ['A', 'B', 'C'];
maxMotifLength = 4;
delta0 = 0.5;
alpha0 = 0.5;

%% Find the optimal model parameters that fits the data the best
[modelParams, optCorr, modelHumanCorr] = optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength, delta0, alpha0);

%% Plot the heatmap
clims = [0.4, 0.9];
imagesc([0.25, 0.95], [0.05, 0.95], modelHumanCorr, clims)
colorbar
xlabel('delta');
ylabel('alpha')