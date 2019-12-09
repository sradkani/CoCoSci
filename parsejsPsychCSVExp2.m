function datatable = parsejsPsychCSVExp2

pilotfiles = dir('Experiment/Results Exp2/exp*');

% initialize csvtable
csvtable = table();

% number of sequences per file
numSeqs = [];

% loop through files and concatenate
for i = 1:length(pilotfiles)
    
    % read file
    currentFile = readtable(pilotfiles(i).name);
    
    numSeqs = [numSeqs sum(diff([0 strcmp(currentFile.test_part, 'test')']) == -1 &...
    ~cellfun(@(x) startsWith(x, '<p>'), currentFile.stimulus'))];
    
    % throw out trial_index and time elapsed (was makign trouble in concatenation)
    currentFile.trial_index = [];
    currentFile.time_elapsed = [];
    
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

lastEventsWithoutExcl = find(diff([0 testIdx']) == -1) - 2;

lastEvents = find(diff([0 testIdx']) == -1 &...
    ~cellfun(@(x) startsWith(x, '<p>'), csvtable.stimulus')) - 2;

[~,~,c] = intersect(lastEvents, lastEventsWithoutExcl);

firstEvents = find(diff([0 testIdx']) == 1);

firstEvents = firstEvents(c);

% indices ranging between events
idx = [];
for i = 1:length(firstEvents)
    idx = [idx firstEvents(i):lastEvents(i)+1];
end


% get rt for key press on last element
datatable.rt = cellfun(@str2num, csvtable(lastEvents, :).rt,'UniformOutput', false);


% get stimulus for each 
stimuli = cellfun(@(x) char(extractBetween(x, '>', '</div>')),...
    csvtable(idx,:).stimulus, 'UniformOutput', false);

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
% get rt for likert scale (filter instruction & non-test stimuli)eq
datatable.RTlikert = cellfun(@str2num,...
    csvtable(lastEvents + 2,:).rt);

% get response to likert scale
datatable.likertscore = cellfun(@(x) extractBetween(x, '"Q0":', '}'),...
    csvtable(lastEvents + 2,:).responses);

%% get info from predicting next item

datatable.RTpredict = cellfun(@str2num,...
    csvtable(lastEvents + 3,:).rt);

datatable.predicted = cellfun(@(x) extractBetween(x, '"Q0":"', '"}'),...
    csvtable(lastEvents + 3,:).responses);

datatable.correct = cellfun(@(x) startsWith(x, 'Correct'),...
    csvtable(lastEvents + 4,:).stimulus);

%% get subject and trial level info

% 2nd argument for repelem: number of sequences for each participant (maybe
% accummulate array when loading in for loop at beginning of function)
datatable.subjID = repelem(unique(csvtable.subjIDs),numSeqs);

% save trial IDs (trials that were stopped such that the same sequences
% appeared, i.e. they couldve gone on to be different)
[~, ~, datatable.trialID] = uniqueRowsCA(datatable.sequences);