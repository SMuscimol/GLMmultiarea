function [PSTH, tRange] = getPSTHfromy(y, dm, trialType, trialDuration, onlycompleted, varargin)
%%% Computes PSTH from data (y) 

if nargin < 5
    onlycompleted = false;
end

if nargin < 4
    trialDuration = 2.;
end

if nargin < 3
    trialType = 'All';
end

% get trial lengths
lengths = dm.dspec.expt.binfun([dm.dspec.expt.trial.duration]);
% get trial indices
if strcmp(trialType,'All')
    trialIdxs = 1:size(lengths,2);
else
    trialIdxs = find(arrayfun(@(j) strcmp(dm.dspec.expt.trial(j).tType,trialType), 1:size([dm.dspec.expt.trial],2)));
end

MAXtrialLength = dm.dspec.expt.binfun(4.);  
tRange = linspace( 0., 4., MAXtrialLength);
if isempty(trialIdxs)
	PSTH = nan(1,MAXtrialLength);
	return
end

PSTH_All = nan(numel(trialIdxs),MAXtrialLength);

% compute PSTH %
tStart = 1;
trialCounter=0;
for i=1:numel(lengths)
    if strcmp(trialType, 'All')
        doTrial=true;
    else
        doTrial = strcmp(trialType , dm.dspec.expt.trial(i).tType);
    end

    if  doTrial
        trialCounter=trialCounter+1;
        PSTH_All (trialCounter,1:lengths(i))=  y(tStart:tStart+lengths(i)-1,1)* 1/dm.dspec.expt.binSize;
    end
    tStart = tStart + lengths(i);
end

% average while removing nans %
PSTH=nanmean(PSTH_All,1);
end