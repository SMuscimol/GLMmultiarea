function [val, cvPart, yTestOut, yP] = fitqpandanalyze(dspec, fitParams, fitOptions, lambda, toRemove, cvPart)

if nargin < 6
	cvPart = [];
end

if nargin < 5
	toRemove = {};
end

% set up the folds %
tTypes = {'Hit','FA','CR','miss','Aborted'};
fractions =  extracttrialfractions(dspec.expt);
allIdxs = [];
Ntrials = numel(dspec.expt.trial);
if isempty(cvPart)
	if fitOptions.allEqual
        % create CV partition objects for each trial type
        for ti=1:numel(tTypes)
            if fractions.(tTypes{ti})~=0
                cvPart.(tTypes{ti}) = cvpartition(int16(fractions.(tTypes{ti})*Ntrials), 'KFold', fitOptions.kFold);
            else
                cvPart.(tTypes{ti}) = [];
            end
        end
    else
        % use a single CV partition object
		cvPart = cvpartition(Ntrials, 'KFold', fitOptions.kFold);
	end
end

for k=1:fitOptions.kFold
    % sample test trials %
    if fitOptions.allEqual
        % sample in a way which balances trial types
        tTest = zeros(Ntrials,1,'logical');
        for ti=1:numel(tTypes)
            if ~isempty(cvPart.(tTypes{ti}))
                allIdxs.(tTypes{ti}) = find(arrayfun(@(tidx) strcmp(dspec.expt.trial(tidx).tType,tTypes{ti}),1:numel(dspec.expt.trial)));
                if k > cvPart.(tTypes{ti}).NumTestSets
                    k2 = randi(cvPart.(tTypes{ti}).NumTestSets);
                    tTest(allIdxs.(tTypes{ti})(test(cvPart.(tTypes{ti}),k2))) = true;
                else
                    tTest(allIdxs.(tTypes{ti})(test(cvPart.(tTypes{ti}), k))) = true;
                end
            end
        end
    else
        % sample using the data with unbalanced fraction of trials
        tTest = test(cvPart, k);
    end

    %%% Define training and testing set. Oversampling is used in training set to balance data %%%
    if fitParams.oversample
        [trainExpt, testExpt, ~] = separateandoversample(dspec.expt, fitParams.finalTestFraction, fractions, tTest);
    end
    
    %%% define training and testing design specifications and matrices %%%
    dspecTrain = dspec;
    dspecTrain.expt = trainExpt;
    dspecTest = dspec;
    dspecTest.expt = testExpt;
    %
    dmTrain = buildGLM.compileSparseDesignMatrix(dspecTrain, 1:size(trainExpt.trial,2));
    dmTest = buildGLM.compileSparseDesignMatrix(dspecTest, 1:size(testExpt.trial,2));
    % z-score the analog vars %
    if fitParams.zscore
        [iS, iE] = getXcovaridxs('whiskSpeed', dmTrain.dspec);
        dmTrain.X(:,iS:iE) = zscore(full(dmTrain.X(:,iS:iE)));
        dmTest.X(:,iS:iE) = zscore(full(dmTest.X(:,iS:iE)));
        [iS, iE] = getXcovaridxs('jawSpeed', dmTrain.dspec);
        dmTrain.X(:,iS:iE) =  zscore(full(dmTrain.X(:,iS:iE)));
        dmTest.X(:,iS:iE) =  zscore(full(dmTest.X(:,iS:iE)));
        [iS, iE] = getXcovaridxs('tongueSpeed', dmTrain.dspec);
        dmTrain.X(:,iS:iE) =  zscore(full(dmTrain.X(:,iS:iE)));
        dmTest.X(:,iS:iE) =  zscore(full(dmTest.X(:,iS:iE)));
    end
    %%% get training and testing target %%%
    yTrain = buildGLM.getBinnedSpikeTrain(trainExpt, 'sptrain', dmTrain.trialIndices);
    yTest = buildGLM.getBinnedSpikeTrain(testExpt, 'sptrain', dmTest.trialIndices);
    yTestOut{k} = yTest;

    % remove regressors if instructed (only for reduced models) %
    [X, ~] = removefromX(dmTrain.X, toRemove, dmTrain.dspec, fitOptions.removalType);

    %%% Fit the model %%%
    % set options %
    opts.family = 'poissexp';
    opts.familyextra = 1;
    opts.getH = true;
    opts.baseline = zeros(length(full(yTrain)),1,class(full(yTrain)));

    %%% set the quadratic form that regularizes and facilitates
    %%% smoothness of regressors
    if strcmp(fitOptions.removalType,'remove')
        qf = qfregressors(dspec, toRemove, fitParams.regrGroups);
    elseif strcmp(fitOptions.removalType,'shuffle')
        qf = qfregressors(dspec, {}, fitParams.regrGroups);
    end
    % add constant column for intercept %
    X = cat(2,ones(size(X,1),1),X);

    %%% run fit %%%
    tmp = glmfitqp(yTrain, X, lambda .* qf, opts);
    
    % store results in the output structure %
    val(k).trainLL = tmp.loglikelihood; % log-likelihood on the training set
    val(k).w = tmp.w; % weights
    val(k).H = tmp.H; % Hessian
    
    %%% test the model %%%
    % remove regressors again, this time from testing set %
    [X, ~] = removefromX(dmTest.X, toRemove, dmTest.dspec, fitOptions.removalType);
    % add column for intercept %
    X = cat(2,ones(size(X,1),1),X);
    % adjust option length %
    opts.baseline = zeros(length(yTest),1,class(yTest));
    
    %%% measure fit quality %%%
    % test Log-Likelihood (LL) %
    val(k).testLL = evalGlmLikelihood(yTest, X, tmp.w, opts.baseline, opts.family, 1.);
    % LL of the constant model %
    totNspikes = sum(full(yTrain)) + sum(full(yTest));
    totTime = ( size(full(yTrain),1) + size(full(yTest),1) ) * fitParams.mainBinsize ; 
    homRate = totNspikes/totTime * fitParams.mainBinsize;
    val(k).homRate = homRate;
    xHom = ones(size(dmTest.X,1),1);
    opts.family = 'poissexp';
    opts.baseline = zeros(length(full(yTest)), 1, class(full(yTest)));
    val(k).homLL = evalGlmLikelihood(full(yTest), xHom, log(homRate), opts.baseline, opts.family, 1.);

    % deviance-based pseudo-R-squared %
    val(k).LLf = evalMaxGlmLikelihood(yTest, 'poissexp', 1);
    val(k).Rsq = (val(k).homLL - val(k).testLL)/(val(k).homLL - val(k).LLf);
    ns = length(yTest); % number of samples
    val(k).adjRsq = 1-(1-val(k).Rsq)*(ns-1)/(ns-1-length(val(k).w)+1); % adjusted R-squared 
    % R-squared for each individual trial %
    val(k).RsqST = zeros(numel(testExpt.trial),1); 
    for t=1:numel(testExpt.trial)
        [iS, iE] = getXtrialidxs(t, testExpt);
        val(k).RsqST(t) = getRfromXy(yTest(iS:iE), X(iS:iE,2:end), val(k).w', yTrain);
    end

    NspTest = sum(full(yTest)); % tot number of spikes in test set
    val(k).MI = log2(exp(1)) .* 1 ./ NspTest .* (val(k).homLL - val(k).testLL); % mutual information regressors/spikes

    % reconstructed PSTH %
    for i=1:numel(tTypes)
        n.wPlain = val(k).w';
        [val(k).PSTH.(tTypes{i}), ~] = getGenPSTH(n, X(:,2:end), testExpt, tTypes{i}, 1., true, false);
        [val(k).PSTHdata.(tTypes{i}),~] = getPSTHfromy(yTest, dmTest, tTypes{i} , 1., true);

    end
    % predicted averaged spike rate %
    yP{k} = exp(X * n.wPlain')./testExpt.binSize;

end

end
