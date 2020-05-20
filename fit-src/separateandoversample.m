function [trainingExpt, testExpt, tTest] = separateandoversample(expt, percentage, fractions, tTest)
%%% Returns an experiment structure for training and one for testing, and
%%% the trial indices used for testing. 
%%% expt : full experiment structure
%%% percentage : fraction of trials used for testing
%%% fractions : fractions of trial for each trial type
%%% tTest : trials indices to use for testing, if externally set

if nargin < 4
	tTest = [];
end

% FIND NUMBER OF TEST TRIALS FOR EACH TYPE %
nTrials = size(expt.trial,2);
nTestTrials = int16(round(percentage * nTrials));
nTest.Hit = int16(floor(nTestTrials * fractions.Hit));
nTest.miss = int16(floor(nTestTrials * fractions.miss));
nTest.CR = int16(floor(nTestTrials * fractions.CR));
nTest.FA = int16(floor(nTestTrials * fractions.FA));
nTest.Aborted = int16(floor(nTestTrials * fractions.Aborted));

% SEPARATE TEST TRIALS (IF NOT PROVIDED) %
types = {'Hit','miss','CR','FA','Aborted'};
if numel(tTest)==0
	for t=1:size(types,2)
    	tmp = find(strcmp({expt.trial.tTypeV},types{t}));
    	tmp = tmp(randperm(size(tmp,2)));
    	toTest = tmp(1:nTest.(types{t}));
    	tTest = horzcat(tTest, toTest);
	end
end

tTrain = 1:nTrials;
tTrain(tTest) = [];

trainingExpt = expt;
testExpt = expt;

trainingExpt.trial(tTest) = [];
testExpt.trial(tTrain) = [];


% FIND MOST REPRESENTED TRIAL TYPE %
nTypes=zeros(1,numel(types));
for t=1:numel(types)
    nTypes(t) = sum(strcmp({trainingExpt.trial.tTypeV},types{t}));
end
nMax=max(nTypes);

% OVERSAMPLING %
for t=1:size(types,2)
    tmp = find(strcmp({trainingExpt.trial.tTypeV},types{t}));
    if (numel(tmp)<nMax) && (numel(tmp)>0)
        toAdd=nMax-numel(tmp);
        trainingExpt.trial = horzcat(trainingExpt.trial,trainingExpt.trial(datasample(tmp,toAdd)));
    end
end

end
