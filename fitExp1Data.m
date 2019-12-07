%% Load and format human data
clear;clc;
parseJsPsychCSVExp1

datatable(strcmp(datatable.stimulus, ''), :) = [];
uniqueIDs = unique(datatable.trialID);
sequences = {};
subjRand_human = {};
for i=1:length(uniqueIDs)  
    seq = datatable.stimulus(datatable.trialID==uniqueIDs(i));
    sequences{i} = seq{1};
    subjRand_human{i} = str2num(datatable.response(datatable.trialID==uniqueIDs(i)));
end



%% Set the model parameters
alphabet = ['A', 'B', 'C'];
maxMotifLength = 4;
delta0 = 0.8;
alpha0 = 0.1;

%% Find the optimal model parameters that fits the data the best
[modelParams, optCorr] = optimizeModelParams(sequences, subjRand_human, alphabet, maxMotifLength, delta0, alpha0);