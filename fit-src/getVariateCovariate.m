function [dm, discTrials] = getVariateCovariate(n, areaIndex, SessionData, fitParams)
%%% This function returns the design matrix structure (design matrix X plus
%%% a the design specification sub-structure) and a list of discarded trial
%%% indices. NOTE: this function is highly dependendent on the experimental
%%% design, and it would need to be completely re-designed for another
%%% experiment/task.

% INITIALIZE THE EXPERIMENT %
expt = buildGLM.initExperiment('s',fitParams.mainBinsize,[],[]);
%%% Define all the necessary variables of the experiment %%%
% Event-like varialbes
expt = buildGLM.registerValue(expt, 'tType','trial Type');
expt = buildGLM.registerValue(expt, 'k', 'trial index in the session system');
expt = buildGLM.registerValue(expt, 'tTypeV', 'trial type according to video data');
expt = buildGLM.registerTiming(expt, 'lickOnsetPre','licking onset preparation');
if fitParams.lickPreLong
    expt = buildGLM.registerTiming(expt, 'lickOnsetPreLong','licking onset preparation - long');
end
expt = buildGLM.registerTiming(expt, 'lickOnsetPost','licking onset');
expt = buildGLM.registerTiming(expt, 'lickOnsetDelay','licking onset');
expt = buildGLM.registerTiming(expt, 'LEDOnset','onset of LED');
expt = buildGLM.registerTiming(expt, 'LEDDelay','delayed LED');
expt = buildGLM.registerTiming(expt, 'soundOnset', 'onset of sound');
expt = buildGLM.registerTiming(expt, 'soundDelay', 'onset of sound');
expt = buildGLM.registerTiming(expt, 'wOnset','whisker onset');
expt = buildGLM.registerTiming(expt, 'wDelay','whisker-related delayed activity');
if ~fitParams.onlyCorrect
    expt = buildGLM.registerTiming(expt, 'reward','reward delivery');
    expt = buildGLM.registerTiming(expt, 'punishmentFA','punishment for false alarms');
end

% continuous variables %
%%% whisking %%%
expt = buildGLM.registerContinuous(expt, 'whiskSpeed','whisking speed');
%%% jaw movements %%%
expt = buildGLM.registerContinuous(expt, 'jawSpeed','Speed of Jaw movements');
%%% tongue movements %%%
expt = buildGLM.registerContinuous(expt, 'tongueSpeed', 'Speed of tongue movements');

% slow variables %
expt = buildGLM.registerContinuous(expt,'tIndex','trial index');
expt = buildGLM.registerContinuous(expt,'prevHit','previous trial was Hit');
expt = buildGLM.registerContinuous(expt,'prevFA','previous trial was FA');
expt = buildGLM.registerContinuous(expt,'prevEL','previous trial was EL');
if fitParams.currentTrialType
    expt = buildGLM.registerContinuous(expt,'tLicking','Licking in current trial');
end

% spike train %
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Neuron spike train');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get the already normalized data %%%
wSpeedNormalized = SessionData.Whisker_sp;
jawSpeedNormalized = SessionData.Jaw_sp;
tongueSpeedNormalized = SessionData.Tongue_sp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTrials = min( size(SessionData.TrialOnsets_All,2) , size(wSpeedNormalized,1) );
kAdded = 1;
N = size(SessionData.Areas(areaIndex).Neurons,2);

durations = arrayfun(@(k) getTrialDuration(k, jawSpeedNormalized, fitParams.videoFrameLength), 1:nTrials);

%%% DATA COLLECTION %%%
discTrials = [];
for k=1:nTrials
    if SessionData.CompletedTrialIndices(k) || fitParams.onlyCompleted==false
        if fitParams.onlyCorrect && (SessionData.HitIndices(k) || SessionData.CRIndices(k)) 
			doTrial = true;
		elseif ~fitParams.onlyCorrect
			doTrial = true;
		else
			doTrial = false;
        end
        
        if doTrial
			duration = durations(k);
        	% If lick is too early --> discard the trial
        	if SessionData.LickOnsets(k) < -1
            	duration = NaN;
        	end
        	if ~isnan(duration)
            	trial = buildGLM.newTrial(expt,duration);
            	
            	% Trial Type %
            	if SessionData.HitIndices(k)
                	trial.tType = 'Hit';
            	elseif SessionData.MissIndices(k)
                	trial.tType = 'miss';
            	elseif SessionData.FAIndices(k)
                	trial.tType = 'FA';
            	elseif SessionData.CRIndices(k)
                	trial.tType = 'CR';
            	else
                	trial.tType = 'Aborted';
            	end
            	%%%%%%%%%%%%%%
            	
            	% Pad or cut video traces and remove NaNs%
            	wSpeedNormalized(isnan(wSpeedNormalized)) = 0.;
            	jawSpeedNormalized(isnan(jawSpeedNormalized)) = 0.;
                tongueSpeedNormalized(isnan(tongueSpeedNormalized)) = 0.;
            	
            	if trial.duration < fitParams.maxTrialDuration
                	jawSpeed = transpose(jawSpeedNormalized(k,1:ceil(trial.duration/fitParams.videoFrameLength)));
                	jawSpeed = cat(1,0.,jawSpeed);
                    %
                	whiskSpeed = transpose(wSpeedNormalized(k,1:ceil(trial.duration/fitParams.videoFrameLength)));
                	whiskSpeed = cat(1,0.,whiskSpeed);
                    %
                    tongueSpeed = transpose(tongueSpeedNormalized(k,1:ceil(trial.duration/fitParams.videoFrameLength)));
                    tongueSpeed = cat(1,0.,tongueSpeed);
                else
                	jawSpeed = transpose(jawSpeedNormalized(k,1:end));
                	jawSpeed = cat(1,0.,jawSpeed);
                    %
                	whiskSpeed = transpose(wSpeedNormalized(k,1:end));
                	whiskSpeed = cat(1,0.,whiskSpeed);
                    %
                    tongueSpeed = transpose(tongueSpeedNormalized(k,1:end));
                    tongueSpeed = cat(1,0.,tongueSpeed);
                end
                
                if fitParams.lickOnsetUseVideo	
            		trial.lickOnsetPre = 2 + SessionData.LickOnsets(k) - fitParams.regressors{fitParams.lickOnsetPreRIndex}.tBefore;
                    if fitParams.lickPreLong
                        trial.lickOnsetPreLong = 2 + SessionData.LickOnsets(k) - fitParams.regressors{fitParams.lickOnsetPreRIndexLong}.tBefore;
                    end
                    trial.lickOnsetPost =  2 + SessionData.LickOnsets(k);
                    trial.lickOnsetDelay = trial.lickOnsetPost + fitParams.delayShift;
                else
					if SessionData.ReactionTimes_All(k) ~= 0 
						trial.lickOnsetPre = 2 + SessionData.ReactionTimes_All(k) - fitParams.regressors{fitParams.lickOnsetPreRIndex}.tBefore;
                        if fitParams.lickPreLong
                            trial.lickOnsetPreLong = 2 + SessionData.ReactionTimes_All(k) - fitParams.regressors{fitParams.lickOnsetPreRIndexLong}.tBefore;
                        end
                        trial.lickOnsetPost = 2 + SessionData.ReactionTimes_All(k);
                        trial.lickOnsetDelay = trial.lickOnsetPost + fitParams.delayShift;
					else
						trial.lickOnsetPre = NaN;
                        if fitParams.lickPreLong
                            trial.lickOnsetPreLong = NaN;
                        end
                        trial.lickOnsetPost = NaN;
                        trial.lickOnsetDelay = NaN;
					end
                end
                
                %%% CONTROLS ON TRIAL DURATION %%%
                if trial.lickOnsetPost > trial.duration
                    disp(strcat('Trial ',num2str(k),' discarded because lick is after end of video'));
                    discTrials = cat(1, discTrials, k);
                    continue
                end
                
                if trial.lickOnsetDelay > trial.duration
                    trial.lickOnsetDelay = NaN;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            	
            	
            	
            	if SessionData.StimIndices_All(k)
                	trial.wOnset = 1. +eps - fitParams.otherShift;
            	else
                	trial.wOnset = NaN;
            	end
            	
            	if isnan(trial.wOnset)
                	trial.wDelay = NaN;
            	else
                	trial.wDelay = min(trial.duration, trial.wOnset + eps + fitParams.regressors{fitParams.wOnsetRIndex}.tAfter);
            	end
            	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            	
            	% Fixed timing events %
            	trial.LEDOnset = 0. + eps;
                trial.LEDDelay = trial.LEDOnset + fitParams.delayShift;
                if trial.duration > 2. + fitParams.delayShift
                	trial.soundOnset = 2. +eps - fitParams.otherShift; % normal versions
                    trial.soundDelay = trial.soundOnset + fitParams.delayShift;
                elseif trial.duration > 2. && (trial.duration <= 2. + fitParams.delayShift)
                    trial.soundOnset = 2. +eps - fitParams.otherShift; % normal versions
                    trial.soundDelay = NaN;
                else
                	trial.soundOnset = NaN;
                    trial.soundDelay = NaN;
                end
            	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            	
            	%%% Slow variables %%%%
            	% Trial index %
            	trial.tIndex = (k / nTrials) .* ones(expt.binfun(trial.duration),1);
            	% Previous trial type
            	trial.prevHit = zeros(expt.binfun(trial.duration),1);
            	trial.prevFA = zeros(expt.binfun(trial.duration),1);
            	trial.prevEL = zeros(expt.binfun(trial.duration),1);
            	if k>1
                	if SessionData.HitIndices(k-1)
                    	trial.prevHit = ones(expt.binfun(trial.duration),1);
                	elseif SessionData.FAIndices(k-1)
                    	trial.prevFA = ones(expt.binfun(trial.duration),1);
                	elseif (~SessionData.MissIndices(k-1)) && (~SessionData.CRIndices(k-1))
                    	trial.prevEL = ones(expt.binfun(trial.duration),1);
                	end
            	end
            	% Current trial type %
            	if fitParams.currentTrialType
                	if ~isnan(trial.lickOnsetPre)
                    	trial.tLicking = ones(expt.binfun(trial.duration),1);
                	else
                    	trial.tLicking = zeros(expt.binfun(trial.duration),1);
                	end
            	end
            	%%%%%%%%%%%%%%%%%%%%%%%%
            	
            	%%% Continuous variables %%%            	
            	% Whisk bouts %
                trial.whiskSpeed = getwhiskerspeed(whiskSpeed, fitParams.binSizeContinuous, trial, expt, fitParams ) ;
            	% jaw movements %
                trial.jawSpeed = getjawspeed(jawSpeed, fitParams.binSizeContinuous, trial, expt, fitParams);                
                % tongue movements %
                trial.tongueSpeed = gettonguespeed(tongueSpeed, fitParams.binSizeContinuous, trial, expt, fitParams);
            	           	
            	%%% reward / sounds delivery %%%
                if ~fitParams.onlyCorrect
                    if strcmp(trial.tType, 'Hit')
                        % check whether there was a valve opening in the trial %
                        onsets = SessionData.ValveOnsets .* fitParams.valveBinsize ...
                            - SessionData.TrialOnsets_All(1,k) ;
                        if ~isempty(find((onsets>=2) .* (onsets<=3.2)))
                            trial.reward = onsets(find((onsets>=2) .* (onsets<=3.2))) - fitParams.otherShift ;
                        else
                            % discard Hit trials without valve opening %
                            disp('Discarded Hit trial because missing valve opening')
                            discTrials = cat(1, discTrials, k);
                            continue
                        end
                        if trial.reward >= trial.duration
                            disp('Hit trial ends before reward!')
                            trial.reward = NaN;
                        end
                    else
                        trial.reward = NaN;
                    end
                    if strcmp(trial.tType, 'FA')
                        trial.punishmentFA = min(trial.duration, max(0.,2. + SessionData.ReactionTimes_All(k) - fitParams.otherShift));
                    else
                        trial.punishmentFA = NaN;
                    end
                end
            	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            	
            	%%% Add trial type based on videos %%%
            	if isnan(trial.wOnset)
                	if isnan(trial.lickOnsetPre)
                    	trial.tTypeV = 'CR';
                	else
                    	if trial.lickOnsetPost > 2
                        	trial.tTypeV = 'FA';
                    	else
                        	trial.tTypeV = 'Aborted';
                    	end
                	end
            	else
                	if isnan(trial.lickOnsetPre)
                    	trial.tTypeV = 'miss';
                	else
                    	if trial.lickOnsetPost > 2
                        	trial.tTypeV = 'Hit';
                    	else
                        	trial.tTypeV = 'Aborted';
                    	end
                	end
                end      	
            	
            	%%% Spike train %%%
            	trial.sptrain = SessionData.Areas(areaIndex).Neurons(n).SpikeTimes(find((SessionData.Areas(areaIndex).Neurons(n).SpikeTimes>SessionData.TrialOnsets_All(k)) & ...
                	SessionData.Areas(areaIndex).Neurons(n).SpikeTimes<SessionData.TrialOnsets_All(k) + trial.duration )) - SessionData.TrialOnsets_All(k);
            	%%%%%%%%%%%%%%%
            	
                % trial index (not discarding removed trials) %
                trial.k = k;
            	expt = buildGLM.addTrial(expt, trial, kAdded);
            	kAdded = kAdded+1;
            else
                if SessionData.HitIndices(k)
                    discTrials = cat(1, discTrials, k);
                end
            end
        end
    end
end

%%% REGRESSOR DESIGN %%%
% design specifications %
dspec = buildGLM.initDesignSpec(expt);

%%% event-like regressors %%%
for i=1:size(fitParams.regressors,2)
    dspec = buildGLM.addCovariateTiming(dspec,fitParams.regressors{i}.covariate, ...
        fitParams.regressors{i}.covariate,fitParams.regressors{i}.description,...
        fitParams.regressors{i}.bs);
end

%%% continuous regressors %%%
% whisker movements
dspec = buildGLM.addCovariateRaw(dspec,'whiskSpeed','effect of whisking speed (bout)');
% jaw movements %
dspec = buildGLM.addCovariateRaw(dspec,'jawSpeed','effect of speed of jaw movements');
% tongue movements %
dspec = buildGLM.addCovariateRaw(dspec, 'tongueSpeed', 'effect of speed of tongue movements');

%%% slow variables %%%
% trial index %
dspec = buildGLM.addCovariateRaw(dspec,'tIndex','effect of time since session start');
% previous trial type %
dspec = buildGLM.addCovariateRaw(dspec,'prevHit','effect of previous trial bein Hit');
dspec = buildGLM.addCovariateRaw(dspec,'prevFA','effect of previous trial bein FA');
dspec = buildGLM.addCovariateRaw(dspec,'prevEL','effect of previous trial bein EL');
% current trial type %
if fitParams.currentTrialType
    dspec = buildGLM.addCovariateRaw(dspec,'tLicking','effect of licking in current trial');
end
%%%%%%%%%%%%%%%%%%%%%%

% Add spike history regressor %
if fitParams.spikeHistory
    bs = basisFactory.makeNonlinearRaisedCos(fitParams.SHnumbins,fitParams.mainBinsize,[0 fitParams.SHeffectlength], 2);
    dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist','sptrain','History Filter',bs);
end

% get the independent variable (design matrix) %
dm = buildGLM.compileSparseDesignMatrix(dspec,1:size(expt.trial,2));

% z-score the analog vars %
if fitParams.zscore
    [iS, iE] = getXcovaridxs('whiskSpeed', dm.dspec);
    dm.X(:,iS:iE) = zscore(full(dm.X(:,iS:iE)));
    [iS, iE] = getXcovaridxs('jawSpeed', dm.dspec);
    dm.X(:,iS:iE) =  zscore(full(dm.X(:,iS:iE)));
    [iS, iE] = getXcovaridxs('tongueSpeed', dm.dspec);
    dm.X(:,iS:iE) =  zscore(full(dm.X(:,iS:iE)));
end


end
