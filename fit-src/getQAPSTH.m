function [PSTH_full, PSTH_Qmod, PSTH_Qdata, PSTH_fullmod] = getQAPSTH(GLMDataArea, results, fitParams, nType)
        
%fitParamsTmp = fitParams;
%fitParamsTmp.zscore = true; % ?

PSTH_Qmod = zeros(size(results.Neurons(1).PSTHgen.Hit));
PSTH_Qdata = zeros(size(results.Neurons(1).PSTHgen.Hit));

ii = 0;
for i=1:numel(results.Neurons)
    doNeuron = true;
    if strcmp(nType,'All')
        doNeuron = true;
    elseif strcmp(nType,'RS')
        doNeuron = results.Neurons(i).PeakToBaseline > 0.34;
    elseif strcmp(nType,'FS')
        doNeuron = results.Neurons(i).PeakToBaseline < 0.26;
    end
    if results.Neurons(i).cvRsq <= 0.01
        doNeuron = false;
    end
    if doNeuron
        ii = ii + 1;
        disp(i)
        n = results.Neurons(i);
        [dm, discTr] = getVariateCovariate(n.n, 1, GLMDataArea.Sessions(n.sIndex), fitParams);

        % divide the session in quiet and active trials
        HitIdxs = find(strcmp({dm.dspec.expt.trial.tType}, 'Hit'));

        %%% check that there is no mismatch in the number of Hit trials %%%
        if numel(find(HitIdxs)) ~= numel(find(GLMDataArea.Sessions(n.sIndex).HitIndices))
            warning('Sessions data and results should have the same number of Hit trials')
        end

        HitIndices = GLMDataArea.Sessions(n.sIndex).HitIndices;
        HitIndices(discTr) = [];

        QAcovar = 'Jaw_sp'; % covariate that determines the split in quiet/active
        QAtrace = GLMDataArea.Sessions(n.sIndex).(QAcovar);
        [Aind , Qind] = FindQA(HitIndices, QAtrace);

        % get model spike counts on quiet trials only 
        yQ = buildGLM.getBinnedSpikeTrain(dm.dspec.expt, 'sptrain', HitIdxs(Qind));
        xB = cat(2,ones(size(dm.X,1),1),full(dm.X));
        yQgen = [];
        for j=HitIdxs(Qind)
            [iS, iE] = getXtrialidxs(j, dm.dspec.expt);
            yQgen = cat(1,yQgen, poissrnd(exp(xB(iS:iE,:) * n.wPlain')));
        end

        % measure model and data PSTHs
        PSTHdata_quiet = nan(numel(find(Qind)),size(n.PSTHdata.Hit,1));
        PSTHgen_quiet = nan(numel(find(Qind)),size(n.PSTHdata.Hit,1));
        try
            Xtmp = mergeXs(n.Xhist, dm);
        catch
            warning('No history regressor present')
            Xtmp = dm.X;
        end
        PSTHcheck = getGenPSTH(n,Xtmp , dm.dspec.expt, 'Hit');
        lengths = dm.dspec.expt.binfun([dm.dspec.expt.trial(HitIdxs(Qind)).duration]);

        tStart = 1;
        trialCounter=0;
        hitCounter = 0;
        for j=1:numel(HitIdxs(Qind))
            trialCounter=trialCounter+1;
            PSTHdata_quiet(trialCounter,1:lengths(j)) =  yQ(tStart:tStart+lengths(j)-1,1)* 1/dm.dspec.expt.binSize;
            PSTHgen_quiet(trialCounter,1:lengths(j))  =  yQgen(tStart:tStart+lengths(j)-1,1)* 1/dm.dspec.expt.binSize;
            tStart = tStart + lengths(j);
        end

        PSTHdata_quiet=nanmean(PSTHdata_quiet,1);
        PSTHgen_quiet=nanmean(PSTHgen_quiet,1);

        PSTH_Qmod(:) = nansum([(ii-1)/ii .* PSTH_Qmod , 1/ii .* PSTHgen_quiet(:)],2);
        PSTH_Qdata(:) = nansum([(ii-1)/ii .* PSTH_Qdata , 1/ii .* PSTHdata_quiet(:)] , 2);
    end
end

PSTH_full = zeros(size(results.Neurons(i).PSTHdata.Hit));
PSTH_fullmod = zeros(size(results.Neurons(i).PSTHgen.Hit));
ii = 0;
for i =1:numel(results.Neurons)
    doNeuron = true;
    if strcmp(nType,'All')
        doNeuron = true;
    elseif strcmp(nType,'RS')
        doNeuron = results.Neurons(i).PeakToBaseline > 0.34;
    elseif strcmp(nType,'FS')
        doNeuron = results.Neurons(i).PeakToBaseline < 0.26;
    end
    if results.Neurons(i).cvRsq <= 0.01
        doNeuron = false;
    end
    if doNeuron
        ii = ii +1;
        disp(i);
        PSTH_full(:)= (ii-1)/ii .* PSTH_full + 1/ii .* results.Neurons(i).PSTHdata.Hit(:);
        PSTH_fullmod(:)= (ii-1)/ii .* PSTH_fullmod + 1/ii .* results.Neurons(i).PSTHgen.Hit(:);
    end
end
end