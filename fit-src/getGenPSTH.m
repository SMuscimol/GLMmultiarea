function [PSTH, tRange] = getGenPSTH(n, X, expt, trialType, trialDuration, onlycompleted, average,varargin)
%%% computes PSTH from a model (neuron object n) and a design matrix %%%

if nargin < 7
	average = false;
end

if nargin < 6
    onlycompleted = true;
end

if nargin < 5
    trialDuration = 1.;
end

if nargin < 4
    trialType = 'All';
end

% concatenate ones for intercept
xB = cat(2,ones(size(X,1),1),full(X));
% sample spikes 
genSpTrain = poissrnd(exp(xB * n.wPlain'));

% compute trial lengths
lengths = expt.binfun([expt.trial.duration]);
% compute trial indices
if strcmp(trialType,'All')
    trialIdxs = 1:size(lengths,2);
else
    trialIdxs = find(arrayfun(@(j) strcmp(expt.trial(j).tType,trialType), 1:size([expt.trial],2)));
end 

MAXtrialLength = expt.binfun(4);  
tRange = linspace( 0., 4., MAXtrialLength);
if isempty(trialIdxs)
	PSTH = nan(1, MAXtrialLength);
	return
end

% compute PSTH %
PSTH_All = nan(numel(trialIdxs),MAXtrialLength);
tStart = 1;
trialCounter = 0;
for i=1:numel(lengths)

    if strcmp(trialType, 'All')
        doTrial=true;
    else
        doTrial = strcmp(trialType , expt.trial(i).tType);
    end

    if doTrial
        trialCounter = trialCounter + 1;   
		if average
			tmp = exp(xB * n.wPlain');
			PSTH_All(trialCounter,1:lengths(i))=  tmp(tStart:tStart+lengths(i)-1,1)* 1/expt.binSize;
		else
        	PSTH_All(trialCounter,1:lengths(i))=  genSpTrain(tStart:tStart+lengths(i)-1,1)* 1/expt.binSize;
		end
    end
    tStart = tStart + lengths(i);
end

PSTH=nanmean(PSTH_All,1);

end
