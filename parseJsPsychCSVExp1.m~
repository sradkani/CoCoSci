pilotfiles = dir('Experiment/Results Exp1/exp*');

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
% get rt for sequences
datatable.RTseq = cellfun(@str2num, csvtable(strcmp(csvtable.test_part, 'test'), :).rt);

% get stimulus for sequence (extract formatting from html string)
datatable.stimulus = cellfun(@(x) extractBetween(x, '>', '</div>'),...
    csvtable(strcmp(csvtable.test_part, 'test'),:).stimulus);

% get stimulus for sequence (remove whitespaces and convert to string)
datatable.stimulus = strrep(datatable.stimulus, ' ', '');

%% get info from likert scale
% get rt for likert scale (filter instruction & non-test stimuli)
datatable.RTlikert = cellfun(@str2num,...
    csvtable(~startsWith(csvtable.stimulus, '<p>') &...
    strcmp(csvtable.test_part, ''),:).rt);

% get response to questionnaire
datatable.response = cell2mat(cellfun(@(x) extractBetween(x, '"Q0":', '}'),...
    csvtable(~startsWith(csvtable.stimulus, '<p>') &...
    strcmp(csvtable.test_part, ''),:).responses));

%% get subject and trial level info
datatable.subjID = csvtable(strcmp(csvtable.test_part, 'test'), :).subjIDs;

[~, ~, datatable.trialID] = uniqueRowsCA(datatable.stimulus);
    
%% save table

%writetable(datatable, 'Experiment/pilotdata.csv')

