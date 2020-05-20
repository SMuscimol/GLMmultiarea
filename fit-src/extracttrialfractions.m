function [f, S] = extracttrialfractions(expt)
% Extracts fractions of neurons for each trial type %

Ntot = numel(expt.trial);
NEL = size(find(strcmp({expt.trial.tType},'Aborted')),2);
NHit = size(find(strcmp({expt.trial.tType},'Hit')),2);
NCR = size(find(strcmp({expt.trial.tType},'CR')),2);
NFA = size(find(strcmp({expt.trial.tType},'FA')),2);
NMiss = size(find(strcmp({expt.trial.tType},'miss')),2);

f.Aborted = NEL/Ntot;
f.Hit = NHit/Ntot;
f.CR = NCR/Ntot;
f.FA = NFA/Ntot;
f.miss = NMiss/Ntot;

S = f.Aborted + f.Hit + f.CR + f.FA + f.miss;

end
