pilotfiles = dir('Experiment/Results Exp2/exp*');

% initialize csvtable
csvtable = table();

% loop through files and concatenate
for i = 1:length(pilotfiles)
    
    % read file
    currentFile = readtable(pilotfiles(i).name);
    
    % add list of subject ids
    currentFile.subjIDs = (repelem(extractBetween(pilotfiles(i).name, 'expdata', '.csv'),...
        size(currentFile, 1)))';
    
    % concatenate tables
    csvtable = [csvtable; currentFile];
    
end

% initialize datatable
datatable = table();

%% get info from test sequences

% all the events where a symbol was displayed 
testIdx = strcmp(csvtable.test_part, 'test');

% get index of last event of each sequence (minus 2 because: -1 appears on
% the first 0 affter a streak, and because we have a blank stimulus after a
% space press which also is registered 
lastEvents = find(diff([0 testIdx']) == -1) - 2;

% get rt for key press on last element
datatable.rt = cellfun(@str2num, csvtable(lastEvents, :).rt,'UniformOutput', false);

% get stimulus for each 
stimuli = cellfun(@(x) char(extractBetween(x, '>', '</div>')),...
    csvtable(testIdx,:).stimulus, 'UniformOutput', false);

% convert into sequences
sequences = cell(1, numel(lastEvents));

% will be index
k = 1;
for i = 1:length(stimuli)
    
      % if stimuli is an empty character it means the sequence has ended and
    % we should move to the next sequence
    if isempty(stimuli{i})
        k = k + 1;
    else
        % concatenate stimuli into sequences
        sequences{k} = [sequences{k} stimuli{i}];
    end
    

    

end

% get stimulus for sequence (remove whitespaces and convert to string)
datatable.sequences = sequences';

%% get info from likert scale
% get rt for likert scale (filter instruction & non-test stimuli)
datatable.RTlikert = cellfun(@str2num,...
    csvtable(strcmp(csvtable.trial_type, 'survey-likert'),:).rt);

% get response to likert scale
datatable.likertscore = cell2mat(cellfun(@(x) extractBetween(x, '"Q0":', '}'),...
    csvtable(strcmp(csvtable.trial_type, 'survey-likert'),:).responses));

%% get info from predicting next item

% TODO 
datatable.RTpredict = cellfun(@str2num,...
    csvtable(strcmp(csvtable.trial_type, '???'),:).rt);

datatable.predicted = cellfun(@str2num,...
    csvtable(strcmp(csvtable.trial_type, '???'),:).rt);





%% get subject and trial level info

% 2nd argument for repelem: number of sequences for each participant (maybe
% accummulate array when loading in for loop at beginning of function)
datatable.subjID = repelem(unique(csvtable.subjIDs), ['???']);

% save trial IDs (trials that were stopped such that the same sequences
% appeared, i.e. they couldve gone on to be different)
[~, ~, datatable.trialID] = uniqueRowsCA(datatable.sequences);