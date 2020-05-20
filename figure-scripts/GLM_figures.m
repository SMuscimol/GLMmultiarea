%%% Samuel P. Muscinelli, 2020 %%%
%% load summarized data
load('GLM.mat'); % assumes that the file is in the same directory

%% definitions
allResultsE = GLM.allResultsE;
allResultsN = GLM.allResultsN;
fitParams = GLM.fitParams;
fitOptions = GLM.fitOptions;

areas = {'wS1', 'wS2', 'wM1', 'A1', 'V1', 'PPC', 'dCA1', 'mPFC', 'Striatum', 'wM2', 'ALM', 'tjM1'};
strExpert = '#800080';
strNovice = '#20b2aa';
colorExpert = sscanf(strExpert(2:end),'%2x%2x%2x',[1 3])/255;
colorNovice = sscanf(strNovice(2:end),'%2x%2x%2x',[1 3])/255;
colorLED = [0.2941    0.6627    0.2000];
colorWhisk = [0.9294    0.5020    0.1098];
colorSound = [0.8863    0.1451    0.1255];
colorLick = [0.1,0.1,0.1];

thRS = 0.34;
thFS = 0.26;

tRange = 0.05:0.1:3.95; % for histogram plots

timeBin = 0.1; % binning of GLM fit

%% discarded neurons
minMI = 0.;
nType = 'RS';
resultsDE.Neurons = filterneurons(allResultsE.Neurons, @(n) n.MI <= minMI);
resultsDN.Neurons = filterneurons(allResultsN.Neurons, @(n) n.MI <= minMI);
resultsDE.allNeurons = allResultsE.Neurons;
resultsDN.allNeurons = allResultsN.Neurons;
if strcmp(nType,'RS')
    resultsDE.Neurons = filterneurons(resultsDE.Neurons, @(n) n.PeakToBaseline>0.34);
    resultsDE.allNeurons = filterneurons(allResultsE.Neurons, @(n) n.PeakToBaseline>0.34);
    resultsDN.Neurons = filterneurons(resultsDN.Neurons, @(n) n.PeakToBaseline>0.34);
    resultsDN.allNeurons = filterneurons(allResultsN.Neurons, @(n) n.PeakToBaseline>0.34);
elseif strcmp(nType,'FS')
    resultsDE.Neurons = filterneurons(resultsDE.Neurons, @(n) n.PeakToBaseline<0.26);
    resultsDE.allNeurons = filterneurons(allResultsE.Neurons, @(n) n.PeakToBaseline<0.26);
    resultsDN.Neurons = filterneurons(resultsDN.Neurons, @(n) n.PeakToBaseline<0.26);
    resultsDN.allNeurons = filterneurons(allResultsN.Neurons, @(n) n.PeakToBaseline<0.26);

end
ls = zeros(length(areas),2);
for a=1:length(areas)
    area = areas{a};
    ls(a,1) = length(filterneurons(resultsDN.Neurons, @(n) strcmp(n.area,area)))/length(filterneurons(resultsDN.allNeurons, @(n) strcmp(n.area,area)));
    ls(a,2) = length(filterneurons(resultsDE.Neurons, @(n) strcmp(n.area,area)))/length(filterneurons(resultsDE.allNeurons, @(n) strcmp(n.area,area)));
    [~, pDiscarded(a)] = prop_test([ls(a,2), ls(a,1)] .* [length(filterneurons(resultsDE.allNeurons, @(n) strcmp(n.area,area))),length(filterneurons(resultsDN.allNeurons, @(n) strcmp(n.area,area)))], ...
            [length(filterneurons(resultsDE.allNeurons, @(n) strcmp(n.area,area))), length(filterneurons(resultsDN.allNeurons, @(n) strcmp(n.area,area)))], false);
end

[~, pDiscardedTot] = prop_test([length(resultsDE.Neurons), length(resultsDN.Neurons)], ...
            [length(resultsDE.allNeurons), length(resultsDN.allNeurons)], false);
b = bar(ls);
b(2).FaceColor = colorExpert;
b(1).FaceColor = colorNovice;
xticklabels(areas)
xtickangle(45)
ylabel({'fraction of', 'discarded neurons'});
set(gca, 'box','off');
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.5, 1.], 'PaperUnits', 'Inches', 'PaperSize', [2.5,1.]);

%% Plot individual regressor modulation across areas (Fig 5C)
pMax = 0.05;
nTypes = {'RS'}; % set to nTypes = {'All', 'RS', 'FS'} to produce plots also for fast-spiking neurons.

filterVar = 'MI';
if strcmp(filterVar,'cvRsq')
    filterValue = 0.01;
elseif strcmp(filterVar,'MI')
    filterValue = 0.;
end

for nT=1:length(nTypes)
    nType = nTypes{nT};

    fractExpert = zeros(length(areas),length(fitOptions.toRemove));
    fractNovice= zeros(length(areas),length(fitOptions.toRemove));
    pFracs =  zeros(length(areas),length(fitOptions.toRemove));

    for a=1:length(areas)
            if strcmp(nType, 'All')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) .* (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            elseif strcmp(nType,'RS')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline > 0.34) * (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline > 0.34) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            elseif strcmp(nType,'FS')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline < 0.26) * (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline < 0.26) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            end
    end

    for r=1:length(fitOptions.toRemove)
        figure('Visible','on');

        % expert and novice
        b = bar(cat(2,fractNovice(:,r),fractExpert(:,r)));
        b(2).FaceColor = colorExpert;
        b(1).FaceColor = colorNovice;

        xticklabels(areas);
        xtickangle(45);
        ylabel({'fraction of', 'modulated neurons'});
        ax = gca;
        ax.FontSize = 8; 
        title(strcat(fitOptions.toRemove{r}{:}));
        
        pvaluefactor = 1.;
        for a=1:numel(areas)
            yLim = ylim();
            if pvaluefactor * pFracs(a,r) < 0.0005
                text(a*1-0.05-0.2,  1. * yLim(2), '***','FontSize',5);
            elseif  pvaluefactor * pFracs(a,r) < 0.005
                text(a*1-0.05-0.1,  1. * yLim(2), '**','FontSize',5);
            elseif  pvaluefactor * pFracs(a,r) < 0.05
                text(a*1-0.05,  1. * yLim(2), '*','FontSize',5);
            end
        end

        set(gca, 'box','off');
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.5, 1.], 'PaperUnits', 'Inches', 'PaperSize', [2.5,1.]);
    end
end

%% Table with all areas/regressors
dx = 0.0175;
dy = 0.04;
pvaluefactor = 1.;

fractExpert = zeros(length(areas),length(fitOptions.toRemove));
fractNovice= zeros(length(areas),length(fitOptions.toRemove));
pFracs =  zeros(length(areas),length(fitOptions.toRemove));

nTypes = {'RS'}; % set to {'RS','FS','All'} to produce tables also for fast-spiking neurons
pMax = 0.05;
minRsq = 0.01;

for nT=1:length(nTypes)
    nType = nTypes{nT};
    %filtering
    filterVar = 'MI'; % set to 'cvRsq' to filter neurons that have CV-R-squared bigger than 0.01
    if strcmp(filterVar,'MI')
        filterValue = 0.;
    elseif strcmp(filterVar,'cvRsq')
        filterValue = 0.01;
    end

    for a=1:length(areas)
            if strcmp(nType, 'All')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) .* (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            elseif strcmp(nType,'RS')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline > 0.34) * (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline > 0.34) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            elseif strcmp(nType,'FS')
                nrnsE = filterneurons(allResultsE.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline < 0.26) * (n.(filterVar) > filterValue) ) );
                fractExpert(a,:) = getfractionmodulated(nrnsE, fitOptions, pMax);
                nrnsN = filterneurons(allResultsN.Neurons, @(n) logical(strcmp(n.area, areas{a}) * (n.PeakToBaseline < 0.26) * (n.(filterVar) > filterValue) ) );
                fractNovice(a,:) = getfractionmodulated(nrnsN, fitOptions, pMax);
                for r=1:length(fitOptions.toRemove)
                    [~, pFracs(a,r)] = prop_test([length(nrnsE), length(nrnsN)] .* [fractExpert(a,r),fractNovice(a,r)], ...
            [length(nrnsE), length(nrnsN)], false);
                end
            end
    end

    figure('Visible','on');
    for a=1:length(areas)
        area = areas{a};
        for r=1:length(fitOptions.toRemove)
            h = subplot(length(areas),length(fitOptions.toRemove),(a-1)*length(fitOptions.toRemove) + r);
            set(h, 'position',[0.1 + (r-1)*dx, 0.95-a*dy, dx, dy]);

            b = bar(1:1:2, [fractNovice(a,r), fractExpert(a,r)], 0.75);
            if pvaluefactor * pFracs(a,r) < 0.0005
                text(1.,  0.9, '***','FontSize',5);
            elseif  pvaluefactor * pFracs(a,r) < 0.005
                text(1.2,  0.9, '**','FontSize',5);
            elseif  pvaluefactor * pFracs(a,r) < 0.05
                text(1.4,  0.9, '*','FontSize',5);
            end
            b.FaceColor = 'flat';
            b.CData(1,:) = colorNovice;
            b.CData(2,:) = colorExpert;
            xticks([]);
            yticks([]);
            ylim([0, 1.]);
            if r==1
                text(-7, 0.3, areas{a})
            end
            if a==length(areas)
                xticks([1.5])
                xticklabels(strcat(fitOptions.toRemove{r}{:}))
                xtickangle(45);
            end
        end
    end
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 14, 8.0], 'PaperUnits', 'Inches', 'PaperSize', [14., 8.0]);
end

%% Venn diagrams (Fig. 5D)
% areas/groups for which the intersection between populations of modulated
% neurons was empty return warning due to the implementation of the Venn
% diagram script (venn.m). However, their correctnedd can be checked manually. 
% Hippocampus/Novice is the only diagram that will not be plotted because
% it has only few neurons sensitive to whisker onset and venn.m does not
% allow to plot only one diagram.

filterVar = 'cvRsq';
if strcmp(filterVar,'cvRsq')
    filterValue = 0.01;
elseif strcmp(filterVar,'MI')
    filterValue = 0.;
end

clc
for a = 1:length(areas)
    
    area = areas{a};
    
    for g=1:2
        for nt={'RS'} % set to {'All','RS','FS'} to produce Venn diagrams for all neuron types
            for normalize=[true]
                if g==1
                    gr = 'Expert';
                    results.Neurons = filterneurons(allResultsE.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)) );
                elseif g==2
                    gr = 'Novice';
                    results.Neurons = filterneurons(allResultsN.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)) );
                end

                if strcmp(nt{1}, 'RS')
                    results.Neurons = filterneurons(results.Neurons, @(n) n.PeakToBaseline > 0.34);
                elseif strcmp(nt{1}, 'FS')
                    results.Neurons = filterneurons(results.Neurons, @(n) n.PeakToBaseline < 0.26);
                end

                for r=1:length(fitOptions.toRemove)
                    allPs.(strcat(fitOptions.toRemove{r}{:})) = arrayfun(@(n) n.pRav(r), results.Neurons);
                end
                
                totN = numel(results.Neurons);
         
                %%% WONSET, WDELAY AND LICKONSETPRE %%%
                try
                    if normalize==false
                        fig = figure('Units','inches','Position',[0,0,4,4]);
                        ax = axes('Units','inches','Position',[0. 0. 4. 4.]);
                    else
                        fig = figure('Units','inches','Position',[0,0,50/sqrt(totN),50/sqrt(totN)]);
                        ax = axes('Units','inches','Position',[0. 0. 50/sqrt(totN), 50/sqrt(totN)]);
                    end
                    axis equal
                    set(gca,'visible','off');
                    
                    xlim([-10,10])
                    ylim([-10,10])
                    xticklabels([])
                    yticklabels([])
                    xticks([])
                    yticks([])              

                   

                    ars = zeros(1,3); 
                    regrs = {'wOnset','wDelay','lickOnsetPre'};
                    for r=1:numel(regrs)
                        ars(r) = numel(find(allPs.(regrs{r})<=0.05));
                    end
                    nzR = find(ars~=0);
                    if numel(nzR)==3
                        intrs = zeros(1,4);
                    elseif numel(nzR)==2
                        intrs = zeros(1);
                    else
                        intrs = [];
                    end
                    r=0;
                    for r1=1:numel(nzR)
                        for r2=(r1+1):numel(nzR)
                            r=r+1;
                            intrs(r) = numel(find((allPs.(regrs{nzR(r1)})<=0.05) .* (allPs.(regrs{nzR(r2)})<=0.05)));
                        end
                    end
                    if numel(nzR)==3
                        intrs(end) = numel(find((allPs.(regrs{nzR(1)})<=0.05) .* (allPs.(regrs{nzR(2)})<=0.05) .* (allPs.(regrs{nzR(3)})<=0.05)));
                    end
                    ars(ars==0) = [];

                    disp(area)
                    disp(strcat('totN:', num2str(totN)))
                    disp(strcat('areas:',num2str(ars)))
                    disp(strcat('inters.:',num2str(intrs)))
                    if numel(ars)==3
                        totArea = sum(ars) - sum(intrs(1:3)) + intrs(4);
                    elseif numel(ars)==2
                        totArea = sum(ars) - intrs(1);
                    elseif numel(ars)==1
                        totArea = ars(1);
                    end
                    [H,S] = myvenn(0.1 .* totArea, ars, intrs);

                    disp(strcat('S areas',num2str(S.CircleArea)))

                catch
                    warning(strcat(area,' - ', gr,' norm:',num2str(normalize),' - ', nt{1}, ' had an error'))
                end
            end
        end
    end
end

%% MI median for all areas - Figure S5E
nTypes = {'RS'}; % set to {'All','RS','FS'} to plot for all neuron types

filterVar = 'MI'; % set to 'cvRsq' to discard neurons whose cv-R-squared is smaller than 0.01
if strcmp(filterVar,'cvRsq')
    filterValue = 0.01;
elseif strcmp(filterVar,'MI')
    filterValue = 0.;
end

for nT=1:numel(nTypes)
    nType = nTypes{nT};

    minRsq = 0.01;
    for a =1:length(areas)
        area = areas{a};
        if strcmp(nType,'All')
            results.Neurons = filterneurons(allResultsE.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue) ) );
            MI.(area) = [results.Neurons.MI]';
            results.Neurons = filterneurons(allResultsN.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)) );
            MINovice.(area) = [results.Neurons.MI]';
        elseif strcmp(nType,'RS')
            results.Neurons = filterneurons(allResultsE.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)*(n.PeakToBaseline>0.34)) );
            MI.(area) = [results.Neurons.MI]';
            results.Neurons = filterneurons(allResultsN.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)*(n.PeakToBaseline>0.34)) );
            MINovice.(area) = [results.Neurons.MI]';
        elseif strcmp(nType, 'FS')
            results.Neurons = filterneurons(allResultsE.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)*(n.PeakToBaseline<0.26)) );
            MI.(area) = [results.Neurons.MI]';
            results.Neurons = filterneurons(allResultsN.Neurons, @(n) boolean(strcmp(n.area,area)*(n.(filterVar)>filterValue)*(n.PeakToBaseline<0.26)) );
            MINovice.(area) = [results.Neurons.MI]';
        end
    end

    tmp = zeros(length(areas),2);
    for a = 1:length(areas)
        tmp(a,1) = median(MINovice.(areas{a}));
        tmp(a,2) = nanmedian(MI.(areas{a}));
    end
    figure('Visible','on')

    % expert and novice %
    b = bar(tmp, 0.75);
    b(1).FaceColor = colorNovice;
    b(2).FaceColor = colorExpert;

    xticklabels(areas);
    xtickangle(45);
    ylabel({'median M.I.', '[bits per spike]'})

    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.5, 1.], 'PaperUnits', 'Inches', 'PaperSize', [2.5, 1.]);
    set(gca, 'box', 'off')
end

%% Quiet / Active model and data histograms (Fig. 5E)
% This analysis can be plotted only for regular-spiking neurons %
gr = 'Expert';
for a=1:numel(areas)
    area = areas{a};
    if strcmp(gr, 'Expert')
        PSTH_full = GLM.QA.(area).PSTH_fullE;
        PSTH_Qmod = GLM.QA.(area).PSTH_QmodE;
        PSTH_Qdata = GLM.QA.(area).PSTH_QdataE;
        PSTH_fullmod = GLM.QA.(area).PSTH_fullmodE;
    elseif strcmp(gr, 'Novice')
        PSTH_full = GLM.QA.(area).PSTH_fullN;
        PSTH_Qmod = GLM.QA.(area).PSTH_QmodN;
        PSTH_Qdata = GLM.QA.(area).PSTH_QdataN;
        PSTH_fullmod = GLM.QA.(area).PSTH_fullmodN;
    end
        
    figure()
    hold all
    plot(tRange, PSTH_fullmod, 'DisplayName','Model-all','LineWidth',2,'Color',[0,0,0])
    plot(tRange, PSTH_full,'LineStyle','--','DisplayName','Data-all','LineWidth',2,'Color',[0.3,0.3,0.3])
    tmp = '#0072bdff';
    c = sscanf(tmp(2:end),'%2x%2x%2x',[1 3])/255; 
    plot(tRange, PSTH_Qmod,'DisplayName', 'Model-quiet','LineWidth',2,'Color',c)
    tmp = '#47b6ffff';
    c = sscanf(tmp(2:end),'%2x%2x%2x',[1 3])/255; 
    plot(tRange, PSTH_Qdata, 'LineStyle','--','DisplayName','Data-quiet','LineWidth',2,'Color',c)
    legend();
    xlabel('t [s]')
    xlim([0,4]);
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.242, 1.129], 'PaperUnits', 'Inches', 'PaperSize', [2.242,1.129]);

end

%% Single neuron panel - PSTHs and betas - Figure S5C %%
s = 10; nn=35; % session / neuron
area = 'ALM';
results.Neurons = filterneurons(allResultsE.Neurons, @(n) strcmp(n.area, area) );
n = results.Neurons(logical(([results.Neurons.sIndex] == s) .* ([results.Neurons.n]==nn) .* ([results.Neurons.n]==nn) ));

figure();
    
subplot(312)

cCRmod = '6868ffff';
cCRmod = sscanf(cCRmod(2:end),'%2x%2x%2x',[1 3])/255;

plot(0.05:0.1:3.95,n.PSTHgen.Hit,'color','black','DisplayName','Hit-model','LineWidth',2)
hold on
plot(0.05:0.1:3.95, n.PSTHdata.Hit,'color','black','LineStyle','--','DisplayName','Hit-data','LineWidth',2)

plot(0.05:0.1:3.95, n.PSTHgen.CR,'color',cCRmod,'DisplayName','CR-model','LineWidth',2)
plot(0.05:0.1:3.95, n.PSTHdata.CR,'color','blue','LineStyle','--','DisplayName','CR-data','LineWidth',2)

xlabel('t [s]');
ylabel('PSTH [Hz]');
xlim([0,4])
set(gca,'Box','off')

subplot(311)
plotweights_exp(n,0.05, 1.,true, true,'LineWidth',2);

ax = gca;
ax.FontSize = 8; 
xlim([0,4])
ylabel('\Delta rate')
set(gca,'Box','off')
% bar inset
axes('Position',[0.8 .8 .075 .075])
%b = bar([(n.psWaldAv.jawSpeed.data<0.05) * n.ws.jawSpeed.data; ...
%  (n.psWaldAv.whiskSpeed.data<0.05) *n.ws.whiskSpeed.data ; ...
%  (n.psWaldAv.tongueSpeed.data<0.05) * n.ws.tongueSpeed.data ]);
m = exp(n.wPlain(1))/timeBin;
b = bar([m*(exp((n.psWaldAv.jawSpeed.data<0.05) * n.ws.jawSpeed.data) - 1), ...
         m*(exp((n.psWaldAv.whiskSpeed.data<0.05) * n.ws.whiskSpeed.data) - 1), ...
         m*(exp((n.psWaldAv.tongueSpeed.data<0.05) * n.ws.tongueSpeed.data) - 1)]);
set(gca, 'box','off')
xlabel('')
ylabel('\beta')

subplot(313)
r = 2;
tType = 'Hit';
P=[]; for i=1:fitOptions.kFold; P = vertcat(P, n.PSTHRgen{r,i}.(tType)); end; P = mean(P,1);
plot(0:0.1:3.9, n.PSTHdata.(tType),'Color','k','LineWidth',2,'LineStyle','--','DisplayName','data')
hold on
plot(0:0.1:3.9,n.PSTHgen.(tType),'LineWidth',1.5, 'Color',0.3 .*[1,1,1], 'DisplayName','full model');
plot(0:0.1:3.9, P,'LineWidth',1.5,'Color',colorWhisk, 'DisplayName','reduced model')
xlabel('t [s]')
ylabel('PSTH [Hz]')
xlim([0.,4.])
set(gca,'Box','off')

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3, 4], 'PaperUnits', 'Inches', 'PaperSize', [3,4]);
%% Single neuron panel - load full data%%

load('GLMfulloutput.mat'); % assumes that the file is in the directory (might take a while), please ensure that you have enough RAM.

%% Trial predictions - Figure S5C

%%% PARAMETER TO BE SET BY THE USER - the defaults correspond to the one in
%%% figure ... ... , but with different semi-randomly chosen trials.
gr = 'Expert'; % mouse group, 'Expert' or 'Novice'
area = 'ALM'; % brain area
s = 10; nn=35; % session / neuron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = GLMfulloutput.(gr).(area);
n = results.Neurons(logical(([results.Neurons.sIndex] == s) .* ([results.Neurons.n]==nn)));
dm = results.Sessions(n.sIndex).dm;

ti = 0;
nTrials = 3;

for t=1:nTrials
    k = t;
    tmp = find(n.RsqST{k} > prctile(n.RsqST{k},50));
    tt = randsample(tmp,1);
    ttExpt = gettexptfromttest(tt, k, n, dm.dspec.expt);
    disp(dm.dspec.expt.trial(ttExpt).tType) % output the trial type
    trs(t) = tt;
    ti = ti +1;
    subplot(1,nTrials, ti)
    plotsingletrialprediction([tt,k], n, dm.dspec.expt,'testSet', 1, 0.1);
    if ti>1
        ylabel('')
    end
end

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 2.756, 1.122], 'PaperUnits', 'Inches', 'PaperSize', [2.756,1.122]);

%% generating heatmap of figure 5C %%
% filter neurons first 
pMax = 0.05;
minMI = 0.;

group = 'Expert';
nType = 'RS';

if strcmp(group, 'Expert')
    results.Neurons = filterneurons(allResultsE.Neurons, @(n) n.MI >= minMI);
elseif strcmp(group, 'Novice')
    results.Neurons = filterneurons(allResultsN.Neurons, @(n) n.MI >= minMI);
end

if strcmp(nType, 'RS')
    results.Neurons = filterneurons(results.Neurons, @(n) n.PeakToBaseline > 0.34);
elseif strcmp(nType, 'FS')
    results.Neurons = filterneurons(results.Neurons, @(n) n.PeakToBaseline < 0.26);
end

% build map %
wNames = fieldnames(results.Neurons(1).ws);
for wName=wNames'
    map.(wName{1}) = zeros(numel(results.Neurons),4);
    map.(wName{1})(:,1) = [results.Neurons.AP];
    map.(wName{1})(:,2) = -[results.Neurons.ML];
    % 3 or 4? leave only the used one %
    map.(wName{1})(:,3) = arrayfun(@(n) mean(exp((n.psWaldAv.(wName{1}).data<pMax) .* n.ws.(wName{1}).data)),  results.Neurons );
    map.(wName{1})(:,4) = arrayfun(@(n) mean((n.psWaldAv.(wName{1}).data<pMax) .* 1/timeBin .* (exp(n.ws.(wName{1}).data + n.wPlain(1))- exp(n.wPlain(1)))),  results.Neurons ); 
end
map.areas = {results.Neurons.area}';
%%% Samuel P. Muscinelli, 2020 %%%
% APPEND CODE TO PLOT MAPS HERE ...% 