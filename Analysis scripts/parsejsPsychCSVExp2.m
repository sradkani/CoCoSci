function datatable = parsejsPsychCSVExp2

alldata = readtable('Experiment2/mTurkExp2/trialdata.csv');

% here the final table will be stored
datatable = table();

% info about each screen shown in jspsych
trialdata = alldata{:, 4}; 

% booleans indicating instructions, trial, prediction question or feedback window
instructions = contains(trialdata, '"test_part": "instructions"');
trials = contains(trialdata, '"test_part": "test"');
prediction = contains(trialdata, '"test_part": "prediction"');
feedback = contains(trialdata, '"test_part": "feedback"');

% get index of last event of each sequence
lastEvents = find(diff([0 trials']) == -1) - 1;

% get index of first event of each sequence
firstEvents = find(diff([0 trials']) == 1);

% indices ranging between events
idx = [];
for i = 1:length(firstEvents)
    idx = [idx firstEvents(i):lastEvents(i)];
end

% get stimulus for each 
stimuli = cellfun(@(x) char(extractBetween(x, '>', '</div>')),...
    trialdata, 'UniformOutput', false);

% get rt for key press on last element
datatable.rt = cellfun(@str2num, csvtable(lastEvents, :).rt,'UniformOutput', false);



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