function datatable = parseJsPsychCSVExp1(excludeOutliers)

alldata = readtable('Experiment1/mTurkExp1/trialdata.csv');

% here the final table will be stored
datatable = table();

% info about each screen shown in jspsych
trialdata = alldata{:, 4}; 

% booleans indicating trial or likert window
botcheck = contains(trialdata, '"test_part": "botcheck"');
trials = contains(trialdata, '"test_part": "test"');
likert = contains(trialdata, '"test_part": "likert"');

% get rt for sequences
datatable.RTseq = cellfun(@(x) extractBetween(x, '"rt": ', ','), trialdata(trials));

% get rt for likert scale (filter instruction & non-test stimuli)
datatable.RTlikert = cellfun(@(x) extractBetween(x, '"rt": ', ','), trialdata(likert));

% get stimulus for sequence (extract formatting from html string)
datatable.stimulus = cellfun(@(x) extractBetween(x, '<div style=\"font-size:40px;\">', '</div>"'), trialdata(trials));

% get stimulus for sequence (remove whitespaces and convert to string)
datatable.stimulus = strrep(datatable.stimulus, ' ', '');

% get response to questionnaire (+1 because 10 is registered as 9 for some
% reason in jspsych)
datatable.response = cellfun(@(x) str2double(extractBetween(x, '"Q0\":', '}')) + 1,...
    trialdata(likert));

% add list of worker ids
[workerIDs, ~, pos] = unique(alldata{:, 1}, 'stable');

trialsPerWorker = pos(trials);

datatable.workerID = workerIDs(trialsPerWorker);

datatable.workerNum = trialsPerWorker;

% check whether participants answered the botcheck right
passedBotCheck = cellfun(@(x) (contains(x, 'green') | contains(x, 'dot')) &...
    (contains(x, 'round') | contains(x, 'circle')), trialdata(botcheck));

% find indices of different workers and fill with botcheck
datatable.passedCheck = passedBotCheck(trialsPerWorker);

[~, ~, datatable.trialID] = uniqueRowsCA(datatable.stimulus);

% mark participants whose average is more than 3 MAD above the median
medianResponse = median(datatable.response);

medianPerWorker = accumarray(datatable.workerNum, datatable.response, [], @median);
MAD = median(abs(medianPerWorker - median(medianPerWorker)));

outlier = medianPerWorker < (medianResponse - 3*MAD) | ...
    medianPerWorker > (medianResponse + 3*MAD);

datatable.outlier = outlier(trialsPerWorker);

if excludeOutliers
    datatable = datatable(~datatable.outlier, :);
    datatable = datatable(~contains(datatable.workerID, 'debug'), :);
end

% keep only those that passed the bot check
datatable = datatable(datatable.passedCheck,:);

save('Exp1_Results.mat', 'datatable')