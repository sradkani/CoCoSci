function datatable =  parsejsPsychCSVExp2

alldata = readtable('Experiment2/mTurkExp2/trialdata.csv');

% here the final table will be stored
datatable = table();

% info about each screen shown in jspsych
trialdata = alldata{:, 4}; 

% booleans indicating instructions, trial, prediction question or feedback window
instructions = contains(trialdata, '"test_part": "instructions"');
trials = contains(trialdata, '"test_part": "test"');
training = contains(trialdata, '"test_part": "training"');
prediction = contains(trialdata, '"test_part": "prediction"');
feedback = contains(trialdata, '"test_part": "feedback"');

% get index of last event of each sequence
lastEvents = find(diff([0 trials']) == -1) - 1;

% get index of first event of each sequence
firstEvents = find(diff([0 trials']) == 1);

% indices ranging between events
sequences = cell(length(firstEvents), 1);
for i = 1:length(firstEvents)
    sequences{i} = char(cellfun(@(x) cell2mat(extractBetween(x, '>', '</div>')),...
    trialdata(firstEvents(i):lastEvents(i),:), 'UniformOutput', false).').';
end

% get stimulus for sequence (remove whitespaces and convert to string)
datatable.sequences = sequences;

% get rt for key press on last element
datatable.rt = cellfun(@(x) cell2mat(extractBetween(x, '"rt": ', ',')),...
    trialdata(lastEvents),'UniformOutput', false);

% get info from predicting next item
datatable.RTpredict = cellfun(@(x) cell2mat(extractBetween(x, '"rt": ', ',')),...
    trialdata(lastEvents + 1), 'UniformOutput', false);

datatable.predicted = cellfun(@(x) extractBetween(x, '{\"Q0\":\"', '\"}"'),...
      trialdata(lastEvents + 1), 'UniformOutput', false);

% check feedback
datatable.correct = cellfun(@(x) contains(x, 'Correct'),...
    trialdata(lastEvents + 2), 'UniformOutput', false);

%% get subject and trial level info

% get worker id's
workerIDs = alldata{:, 1};

% 2nd argument for repelem: number of sequences for each participant (maybe
% accummulate array when loading in for loop at beginning of function)
datatable.workerID = workerIDs(lastEvents);

% save trial IDs (trials that were stopped such that the same sequences
% appeared, i.e. they couldve gone on to be different)
[~, ~, datatable.trialID] = uniqueRowsCA(datatable.sequences);