function datatable = parseJsPsychCSVExp1

alldata = readtable('Experiment1/mTurkExp1/trialdata.csv');

% initialize csvtable
subjTable = table();
    
% add list of worker ids
subjTable.workerIDs = unique(alldata{:, 1});

% here the final table will be stored
datatable = table();

% info about each screen shown in jspsych
trialdata = alldata{:, 4}; 

% booleans indicating instructions, trial or likert window
instructions = contains(trialdata, '"trial_type": "text"') | contains(trialdata, '<p>Thanks for participating');
trials = contains(trialdata, '"stimulus": "<div style=');
likert = contains(trialdata, '"trial_type": "survey-likert"');

% get rt for sequences
datatable.RTseq = cellfun(@(x) extractBetween(x, '"rt": ', ','), trialdata(trials));

% get rt for likert scale (filter instruction & non-test stimuli)
datatable.RTlikert = cellfun(@(x) extractBetween(x, '"rt": ', ','), trialdata(likert));

% get stimulus for sequence (extract formatting from html string)
datatable.stimulus = cellfun(@(x) extractBetween(x, '<div style=\"font-size:40px;\">', '</div>"'), trialdata(trials));

% get stimulus for sequence (remove whitespaces and convert to string)
datatable.stimulus = strrep(datatable.stimulus, ' ', '');

%% get info from likert scale
% get rt for likert scale (filter instruction & non-test stimuli)
datatable.RTlikert = cellfun(@str2num,...
    csvtable(~startsWith(csvtable.stimulus, '<p>') &...
    strcmp(csvtable.test_part, ''),:).rt);

% get response to questionnaire
datatable.response = cellfun(@(x) extractBetween(x, '"Q0":', '}'),...
    csvtable(~startsWith(csvtable.stimulus, '<p>') &...
    strcmp(csvtable.test_part, ''),:).responses);

% insert nan
datatable.response(cellfun(@(x) strcmp(x,'""'), datatable.response)) = {'nan'};


%% get subject and trial level info
datatable.subjID = csvtable(strcmp(csvtable.test_part, 'test'), :).subjIDs;

[~, ~, datatable.trialID] = uniqueRowsCA(datatable.stimulus);
    
%% save table

%writetable(datatable, 'Experiment/pilotdata.csv')

