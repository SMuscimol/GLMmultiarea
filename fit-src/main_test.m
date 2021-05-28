%%% Samuel Pavio Muscinelli, 2020 %%%
%%% This file fits a GLM encoding model to each neuron in the dataset. It
%%% also fits several reduced models to run LR tests. All this
%%% is done for both Novice mice (MICEGROUP=0) and expert mice
%%% (MICEGROUP=1).
%%% The <areas> variable should be set BEFORE running this script, and can contain one or more brain areas for which the fit has to be run.
%%% Considering that it takes approximately one minute to fit each neuron, and that there are on the order of hundred neuron per area, fitting one area can take several hours.
%%% After having fit all the neurons, you can run the script "collect_all_results.m" to reformat the results in a single file like the one used to generate figures.
%%% WARNING %%%
%%% Occasionally (only <10 times for the whole fit), the evidence
%%% optimization algorithm gets caught in a loop, caused by problem with
%%% the package fitglmqp that I am relying on and that I did not try to
%%% solve. If this happens, the easiest way out is to stop the fit, and
%%% then setup the loop so to skip that neurons. In the current version, only two neurons in the whole dataset were skipped.

%%% Set constants %%%
WORK_HOME = getenv('WORK_HOME');  % WARNING: the environment variable WORK_HOME has to be set (you can do that using the function setenv)
testing = false;  % set to true to fit only neuron and verify the whole pipeline.
if testing
    areas = {'wS1'};
else
    if exist('areas', 'var') == 0  % if the area was not provided
	areas = {'wS1', 'wS2', 'wM1', 'A1', 'V1', 'PPC', 'dCA1', 'mPFC', 'Striatum', 'wM2', 'ALM', 'tjM1'};
    end
end

setparamsfile = 'set_parameters.m'; % file that sets parameters
setregressorsfile = 'set_regressors.m'; % file that defines the regressors
gitDir = strcat(WORK_HOME, '/GLMmultiarea'); % directory of the git repo
dataDir = '/burg/theory/users/spm2176/GDplusMemory/'; % directory where the data is
dataName = 'GLMDataV9'; % string (name of the data file, without .mat and without the area-specific part.
outName = 'out'; % name of the output file, without '.mat'

% create output directory if needed %
outDir = strcat(dataDir,date,filesep);
try
		mkdir(outDir);
catch
end

%%% Import modules %%%
addpath(strcat(gitDir,filesep,'neuroGLM',filesep));
addpath(strcat(gitDir,filesep,'fitglmqp',filesep));

if testing
    micegroups = [1];
else
    if exist('micegroups', 'var') == 0 % if the group was not provided
	micegroups = [0, 1];
    end
end

for MICEGROUP = micegroups
    % init filename array %
    if MICEGROUP==1
        files = cell(0);
    elseif MICEGROUP==0
        filesNovice = cell(0);
    end

    for a=1:length(areas)

        area = areas{a};
        disp(area)
        % remove unnecessary variables to empty memory and avoid conflicts %
        clearvars -except testing MICEGROUP micegroups area areas setparamsfile setregressorsfile gitDir dataDir dataName outName outDir figDir nameRoot files filesNovice

        %%% LOAD DATA %%%
        load(strcat(dataDir,dataName,area,'.mat'),'GLMDataArea');
        GLMData = GLMDataArea;
        clear GLMDataArea;
        %%%%%%%%%%%%%%%%%%

        %%% SET RANDOM SEED %%%
        t = datevec(now);
        rng('shuffle');

        %%% load all the parameters %%%
        run(setparamsfile);

        %%% initialize experiment just to construct the basis functions %%%
        expt = buildGLM.initExperiment('s',fitParams.mainBinsize,[],[]);

        %%% load the regressors %%%
        run(setregressorsfile);

        %%%%%%%%%%%%%%%%%%%
        %%% Define function that filters sessions based on mice group %%%%
        filtersession = @(SessionData) SessionData.group == MICEGROUP;
        %%%%%%%%%%%%%%%%%%%

        % name of the output data files %
        toSave = strcat(outDir,outName,area,num2str(MICEGROUP),'_',num2str(rand()),'.mat');  % main output
        toSaveQA = strcat(outDir, outName, 'QA', area, num2str(MICEGROUP),'_',num2str(rand()),'.mat');  % output of Quiet / Active analysis

        import   buildGLM.*

        nSessions = numel(GLMData.Sessions);

        %%% Initialize output struct %%%
        results.datafile = strcat(dataDir,dataName,area,'.mat');
        results.area = area;
        results.Neurons = [];

        % loop over sessions %
        tmp = find([GLMData.Sessions.group]==MICEGROUP);
        if testing
	    maxSess = 1;
        else
	    maxSess = numel(tmp);
        end
	for sIndex=tmp(1:maxSess)
	    disp(sIndex)
            if (filtersession(GLMData.Sessions(sIndex))) % filter session depending on mice group
                
                nMax = size(GLMData.Sessions(sIndex).Areas(1).cluster,2);

                dm=[];
                y = [];

                i=0;
                tTest = [];

                % loop over neurons %
		if testing
		    nMax = 1;
	        end
                for n=1:nMax
		    if (~(strcmp(area, 'V1') & sIndex == 8 & n == 10)) & (~(strcmp(area, 'mPFC') & sIndex == 22 & n == 52))  % exlude two neurons for which EO fails.
                    i=i+1;             
                    disp([sIndex,n])
                    %%% get independent variables from input data
                    dm = getVariateCovariate(n, 1, GLMData.Sessions(sIndex), fitParams);

                    % Save the spike-history part of X %
                    if ~isempty(find(strcmp({dm.dspec.covar.label},'hist')))
                        [~, Xhist] = removefromX(dm.X, {'hist'}, dm.dspec);
                    else
                        Xhist = [];
                    end
                    Neuron.Xhist = Xhist;

                    %%% perform fit %%%
                    % check for NaNs
                    if isempty(find(isnan(dm.X)))

                        %%% perform evidence optimization (EO) to choose lambda %%%
                        y = buildGLM.getBinnedSpikeTrain(dm.dspec.expt, 'sptrain', dm.trialIndices);

                        qfsEO = qfregressors(dm.dspec, {}, fitParams.regrGroups);
                        XEO = cat(2,ones(size(full(dm.X),1),1),full(dm.X));
                        % options for EO %
                        optsEO.family = 'poissexp';
                        optsEO.familyextra = 1;
                        optsEO.getH = true;
                        optsEO.minlambda = fitOptions.lambdaMin;
                        optsEO.maxlambda = fitOptions.lambdaMax;
                        optsEO.lambdaeps = 1e-4;
                        optsEO.lambda0 = 1e-4;
                        EO = evidenceglmfitqp(full(y),XEO,qfsEO,optsEO); % run EO
                        lambda = EO.lambda; % lambda that it is going to be used for the fit

                        %%% train full model %%%
                        [val, cvPart, yTest, yP] = fitqpandanalyze(dm.dspec, fitParams, fitOptions, lambda);				

                        %%% averaging across folds %%%
                        % prepare arrays for results %
                        tmpW = zeros(fitOptions.kFold, size(val(1).w(1:end),1));
                        RsqST = cell(fitOptions.kFold,1);
                        pWald = zeros(fitOptions.kFold,length(val(1).w));
                        Hs = cell(fitOptions.kFold,1);
                        statWald  =  zeros(fitOptions.kFold, length(val(1).w));
                        for mInd=1:fitOptions.kFold
                            tmpW(mInd,:) = val(mInd).w(1:end); % model weights for a single fold
                            RsqST{mInd} = val(mInd).RsqST; % single-trial pseud-R-squared
                            % Wald test %
                            nTmp.wPlain = tmpW(mInd,:); % model weights
                            nTmp.H = val(mInd).H; % model Hessian
                            Hs{mInd} = nTmp.H;
                            [pWald(mInd,:), statWald(mInd,:)] = mywaldtest(nTmp, dm); % returns p-values and relevant stat for Wald test
                        end
                        avW = mean(tmpW,1); % weights averaged across folds
                        cvRsq = mean([val(:).Rsq]); % cross-validated R-squared
                        cvAdjRsq = mean([val(:).adjRsq]); % adjustd cross-validated R-squared
                        MI = mean([val(:).MI]); % mutual information between input and spikes
                        pWaldAv = 1 - chi2cdf( mean(statWald,1), 1); % p-value obtained from Wald test after averaging the relevant statistics
                        psWaldAv = combinePs(dm, pWaldAv(2:end), avW(2:end)); % as above, but recombined to be easily inspected

                        % average test PSTHs
                        tTypes = {'Hit','FA','CR','miss','Aborted'}; % trial types.
                        for ti = 1:numel(tTypes)
                            % model and data PSTHs
                            PSTHgen.(tTypes{ti}) = nan( min(arrayfun(@(kf) numel(val(kf).PSTH.(tTypes{ti})),1:numel(val))),1 );
                            PSTHdata.(tTypes{ti}) = nan( min(arrayfun(@(kf) numel(val(kf).PSTHdata.(tTypes{ti})),1:numel(val))),1 );
                            for kf = 1:fitOptions.kFold
                                if ~isempty(val(kf).PSTH.(tTypes{ti}))
                                    PSTHgen.(tTypes{ti}) =  nansum([(kf-1)/kf .* PSTHgen.(tTypes{ti})  1/kf .* val(kf).PSTH.(tTypes{ti})(1:numel(PSTHgen.(tTypes{ti})))'],2 );
                                    PSTHdata.(tTypes{ti}) = nansum([(kf-1)/kf .* PSTHdata.(tTypes{ti})  1/kf .* val(kf).PSTHdata.(tTypes{ti})(1:numel(PSTHdata.(tTypes{ti})))'] ,2);
                                end
                            end
                        end

                        %%% train reduced models -- list procedure %%%
                        % define result structures - they are largely the same as for the full model %
                        avWR = cell(numel(fitOptions.toRemove),1);
                        cvRsqR = cell(numel(fitOptions.toRemove),1);
                        cvAdjRsqR = cell(numel(fitOptions.toRemove),1);
                        cvPartR = cell(numel(fitOptions.toRemove),1);
                        MIR = cell(numel(fitOptions.toRemove),1);
                        cvLLR =  zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        trainLLR =  zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        LR = zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        LRtrain = zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        pR = zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        pRtrain = zeros(numel(fitOptions.toRemove),fitOptions.kFold);
                        kR = zeros(length(fitOptions.toRemove));
                        PSTHRgen = cell(numel(fitOptions.toRemove),fitOptions.kFold);


                        if fitOptions.removesome
                            % loop over reduced models %
                            for r=1:length(fitOptions.toRemove)
                                % fit reduced model %
                                [valR, cvPartR{r},~,~] = fitqpandanalyze(dm.dspec, fitParams, fitOptions, lambda, fitOptions.toRemove{r}, cvPart);

                                % averaging across folds %
                                tmpWR = zeros(fitOptions.kFold, size(valR(1).w(1:end),1));
                                for mInd=1:fitOptions.kFold
                                    tmpWR(mInd,:) = valR(mInd).w(1:end);
                                    cvLLR(r,mInd) = valR(mInd).testLL;
                                    trainLLR(r,mInd) = valR(mInd).trainLL; 
                                    LR(r,mInd) = 2*(-val(mInd).testLL + valR(mInd).testLL); % likelihood ratio
                                    LRtrain(r,mInd) = 2*(-val(mInd).trainLL + valR(mInd).trainLL); % likelihood ratio on training set
                                    kR(r) = length(val(mInd).w) - length(valR(mInd).w); % difference in number of regressors
                                    pR(r,mInd) = 1 - chi2cdf(LR(r,mInd),kR(r)); % p-value from LR test
                                    pRtrain(r,mInd) = 1 - chi2cdf(LRtrain(r,mInd), kR(r)); % p-value from LR test, on training set
                                end
                                avWR{r} = mean(tmpWR,1);
                                cvRsqR{r} = mean([valR(:).Rsq]);
                                cvAdjRsqR{r} = mean([valR(:).adjRsq]);
                                MIR{r} = mean([valR(:).MI]);
                                avLR(r) = mean(LR(r,:));
                                pRav(r) = 1 - chi2cdf(avLR(r), kR(r));
                                pRtrainAv(r) = 1 - chi2cdf(mean(LRtrain(r,:)), kR(r));
                                % average test PSTHs
                                tTypes = {'Hit','FA','CR','miss','Aborted'};
                                for ti = 1:numel(tTypes)                        
                                    % average across folds
                                    for kf = 1:fitOptions.kFold
                                        if ~isempty(val(kf).PSTH.(tTypes{ti}))
                                            % keep all folds
                                            PSTHRgen{r,kf}.(tTypes{ti}) = valR(kf).PSTH.(tTypes{ti});
                                        end
                                    end
                                end
                            end
                        end

                        %%% Store in the final result structure %%%
                        Neuron.wPlain = avW(1:end);
                        Neuron.ws = buildGLM.combineWeights(dm, avW(2:end));
                        Neuron.cvLL = [val.testLL];
                        Neuron.trainLL = [val.trainLL];
                        Neuron.cvLLR = cvLLR;
                        Neuron.trainLLR = trainLLR;
                        Neuron.homLL = [val.homLL];
                        Neuron.lambda = lambda;
                        Neuron.cvRsq = cvRsq;
                        Neuron.cvAdjRsq = cvAdjRsq;
                        Neuron.cvRsqR = cvRsqR;
                        Neuron.cvAdjRsqR = cvAdjRsqR;
                        Neuron.RsqST = RsqST;
                        Neuron.MI = MI;
                        Neuron.MIR = MIR;
                        Neuron.LR = LR;
                        Neuron.LRtrain = LRtrain;
                        Neuron.avLR = avLR;
                        Neuron.pR = pR;
                        Neuron.pRtrain = pRtrain;
                        Neuron.pRav = pRav;
                        Neuron.pRtrainAv = pRtrainAv;
                        % Wald results
                        Neuron.H = Hs;
                        Neuron.pWaldAv = pWaldAv;
                        Neuron.psWaldAv = psWaldAv;
                        Neuron.pWald = pWald;
                        Neuron.statWald = statWald;
                        %%% test
                        Neuron.PSTHRgen = PSTHRgen;
                        Neuron.tmpW = tmpW;
                        Neuron.PSTHgen = PSTHgen;
                        Neuron.PSTHdata = PSTHdata;
                        % quantities useful to reconstruct spike trains / PSTHs %
                        Neuron.cvPart = cvPart;
                        Neuron.cvPartR = cvPartR;
                        Neuron.yTest = yTest;
                        Neuron.yP = yP;

                        try
                            Neuron.avWR = avWR;
                        catch
                            warning('no reduced model present!')
                        end

                        % saving session-cluster info %
                        Neuron.cluster = GLMData.Sessions(sIndex).Areas(1).cluster(n);
                        Neuron.n= n;
                        Neuron.mouse = GLMData.Sessions(sIndex).mouse;
                        Neuron.date = GLMData.Sessions(sIndex).date;
                        Neuron.group = GLMData.Sessions(sIndex).group;
                        Neuron.sess = GLMData.Sessions(sIndex).sess;
                        Neuron.sIndex = sIndex;
                        Neuron.ClusterDepth = GLMData.Sessions(sIndex).Areas(1).Neurons(n).ClusterDepth;
                        Neuron.PeakToBaseline = GLMData.Sessions(sIndex).Areas(1).Neurons(n).PeakToBaseline;
                        % store additional neuron information, if present %
                        try
                            Neuron.MDS = GLMData.Sessions(sIndex).MDS;
                            Neuron.ARAindex = GLMData.Sessions(sIndex).Areas(1).Neurons(n).ARAindex;
                            Neuron.struct = GLMData.Sessions(sIndex).Areas(1).Neurons(n).struct;
                            Neuron.struct_acr = GLMData.Sessions(sIndex).Areas(1).Neurons(n).struct_acr;
                            Neuron.AP = GLMData.Sessions(sIndex).Areas(1).Neurons(n).AP;
                            Neuron.ML = GLMData.Sessions(sIndex).Areas(1).Neurons(n).ML;
                        catch
                        end

                        results.Neurons = cat(1, results.Neurons, Neuron);

                    else
                        disp('There were some NaNs left in dm.X - this neuron will be discarded')
                    end
	    	end

                end
                results.Sessions(sIndex).dm = dm;
                % save session data %
                save(toSave,'results','fitOptions','fitParams','-v7.3');
            end
        end
        % update filename array %
        if MICEGROUP==1
            [PSTH_full, PSTH_Qmod, PSTH_Qdata, PSTH_fullmod] = getQAPSTH(GLMData, results, fitParams, 'RS');
            QA.(area).PSTH_fullE = PSTH_full;
            QA.(area).PSTH_QmodE = PSTH_Qmod;
            QA.(area).PSTH_QdataE = PSTH_Qdata;
            QA.(area).PSTH_fullmodE = PSTH_fullmod;
            files = cat(2, files, toSave);
        elseif MICEGROUP==0
            [PSTH_full, PSTH_Qmod, PSTH_Qdata, PSTH_fullmod] = getQAPSTH(GLMData, results, fitParams, 'RS');
            QA.(area).PSTH_fullN = PSTH_full;
            QA.(area).PSTH_QmodN = PSTH_Qmod;
            QA.(area).PSTH_QdataN = PSTH_Qdata;
            QA.(area).PSTH_fullmodN = PSTH_fullmod;
            files = cat(2, filesNovice, toSave);
        end
	save(toSaveQA, 'QA', 'files', '-v7.3');
    end
end

%%% Samuel Pavio Muscinelli, 2020 %%%
