function datatable =  parsejsPsychCSVExp2

alldata = readtable('Experiment2/mTurkExp2/trialdata.csv');

% here the final table will be stored
datatable = table();

% info about each screen shown in jspsych
trialdata = alldata{:, 4}; 

% booleans indicating instructions, trial, prediction question or feedback window
botcheck = contains(trialdata, '"test_part": "botcheck"');
trials = contains(trialdata, '"test_part": "training"');
betweenseqs = contains(trialdata, '"test_part": "between_sequences"');

% get index of last event of each sequence
lastEvents = find(diff([0 trials']) == -1) - 1;

% get index of first event of each sequence
firstEvents = find(diff([0 trials']) == 1);

% indices ranging between events
sequences = cell(length(firstEvents), 1);
for i = 1:length(firstEvents)
    sequences{i} = cellfun(@(x) cell2mat(extractBetween(x, 'Slide', '.png')),...
    trialdata(firstEvents(i):lastEvents(i),:), 'UniformOutput', false).';
end

alphabet = 'ABC';
convertedSeqs = cell(size(sequences));
% convert into ABC for further analysis
for i = 1:length(sequences)
    [~,~, uniquevals] = unique(sequences{i,:});
    convertedSeqs{i} = char(alphabet(uniquevals));
end
 


% get stimulus for sequence (remove whitespaces and convert to string)
datatable.sequences = convertedSeqs;

% get rt for key press on last element
datatable.rt = cellfun(@(x) cell2mat(extractBetween(x, '"rt": ', ',')),...
    trialdata(lastEvents),'UniformOutput', false);

%% get subject and trial level info

% get worker id's
workerIDs = alldata{:, 1};

% insert workerID
datatable.workerID = workerIDs(lastEvents);

% check whether participants answered the botcheck right
passedBotCheck = cellfun(@(x) (contains(x, 'green') | contains(x, 'dot')) &...
    (contains(x, 'round') | contains(x, 'circle')), trialdata(botcheck));

% find indices of different workers and fill with botcheck
%passedCheck = repelem(passedBotCheck, diff([find(botcheck), length(botcheck) + 1])).';

% datatable.passedCheck = passedCheck(lastEvents);

% save trial IDs (trials that were stopped such that the same sequences
% appeared, i.e. they couldve gone on to be different)
[~, ~, datatable.trialID] = uniqueRowsCA(datatable.sequences);

% keep only those that passed the check
% datatable = datatable(datatable.passedCheck,:);

datatable = datatable(cellfun(@length, datatable.sequences) < 28, :);
datatable = datatable(cellfun(@length, datatable.sequences) > 3, :);


save('Exp2_Results.mat', 'datatable')