%%% OPTIONS %%%
fitOptions.crossValidated = true;
fitOptions.kFold = 5;
fitOptions.regType = 'poisson';
fitOptions.ridge = true;
fitOptions.lambda = 0.01;
fitOptions.lambdaMin = 1e-5; 
fitOptions.lambdaMax = 1e10;
fitOptions.removesome = true; % whether remove some regressors
fitOptions.removalType = 'remove'; % how to 'remove' some regressors. Can be 'remove' or 'shuffle'
fitOptions.allEqual = true; % whether to maintain the fractions of trial types in the test sets
% regressor or set of regressors to be removed in reduced models %
fitOptions.toRemove{1} = {'wOnset'};
fitOptions.toRemove{2} = {'wDelay'};
fitOptions.toRemove{3} = {'lickOnsetPre'};
fitOptions.toRemove{4} = {'wDelay','lickOnsetPre'};
fitOptions.toRemove{5} = {'jawSpeed'};
fitOptions.toRemove{6} = {'whiskSpeed'};
fitOptions.toRemove{7} = {'tongueSpeed'};
fitOptions.toRemove{8} = {'tongueSpeed', 'jawSpeed'};
fitOptions.toRemove{9} = {'jawSpeed','whiskSpeed','tongueSpeed'};
fitOptions.toRemove{10} = {'LEDOnset'};
fitOptions.toRemove{11} = {'LEDDelay'};
fitOptions.toRemove{12} = {'soundOnset'};
fitOptions.toRemove{13} = {'soundDelay'};
fitOptions.toRemove{14} = {'lickOnsetPost'};
fitOptions.toRemove{15} = {'lickOnsetDelay'};
fitOptions.toRemove{16} = {'tIndex'};
fitOptions.toRemove{17} = {'prevHit'};
fitOptions.toRemove{18} = {'prevFA'};
fitOptions.toRemove{19} = {'prevEL'};
%%%%%%%%%%%%%%

%%% FIT OPTIONS %%%
if strcmp(fitOptions.regType,'poisson')
    fitParams.mainBinsize = 0.1; % seconds ; coarser binning for testing
elseif strcmp(fitOptions.regType,'binomial')
    fitParams.mainBinsize = 0.002; % test for logistic reg.
end
fitParams.maxTrialDuration = 4.0; % seconds
fitParams.onlyCompleted = false;
fitParams.onlyCorrect = true;
fitParams.lickPreLong = false;
fitParams.lickOnsetUseVideo = false; % use video data to detect lick onset
fitParams.lickVideoFrom = 'tongue'; % can be 'jaw' or 'tongue' depending on which video we want to use to detect lick onset
fitParams.tongueThreshold = 0.9;
fitParams.jawThreshold = 0.2;
fitParams.whiskThreshold = 0.008;
fitParams.videoFrameLength = 0.002;
fitParams.videoFrameMultiplier = 10.;
fitParams.spikeHistory = false;
fitParams.SHnumbins = 5; % number of bins in the spike history regressor
fitParams.SHeffectlength = 0.3; % seconds. Duration of spike history effect.
fitParams.currentTrialType = false;
fitParams.valveBinsize = 1.; % seconds - for newer data
fitParams.finalTestFraction = 0.2; % fraction of trials that go in the final test set
fitParams.oversample = true; % choose whether to duplicate samples to balance trial types.
fitParams.fixedtestset = true; %use the same test set across one session
% z-score analog variables %
fitParams.zscore = true;

%%%% regressor shifts backwards %%%%
fitParams.lickBaseRegressorShift = 1.; % backshift of the lick onset regressor
fitParams.lickBoutRegressorShift = 0.1;
fitParams.delayShift = 0.2; % forward shift of delay regressors
fitParams.otherShift = 0.; % basckshift of most of the other regressors
fitParams.binSize = 0.1; % seconds - this is the size of a regressor bin
if strcmp(fitOptions.regType,'poisson')
    fitParams.binSizeContinuous = 0.1; % binSize for the continuous regressors
elseif strcmp(fitOptions.regType,'binomial')
    fitParams.binSizeContinuous = 0.002; % binSize for the continuous regressors
end
%%%
fitParams.regressors = cell(1,1);
%
mainBoxSize = 0.1; %seconds

