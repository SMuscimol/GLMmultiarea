%% collect all results in a lightweight structure 'GLM' %%
DATADIR = 'DIR_WHERE_GLMFIT_OUTPUT_IS';  % ADJUST TO YOUR NEED

% 'VERSION' can be set to either 'light' or 'heavy'.
% The light version produces a structure analogous to the one of GLM.mat, which can be used to generate most figure, with the exception of single trial predictions.
% 'heavy' instead produces an output which is roughly double the size and that also contains single trial information. It is the analogous of the GLMfulloutput.mat file.
VERSION = 'heavy';

tmp = dir(DATADIR);
tmp = tmp(arrayfun(@(i) contains(tmp(i).name, 'outQA'), 1:numel(tmp))); % remove stuff that is not the out files
tmpNovice = tmp(arrayfun(@(i) contains(tmp(i).name, '0_'), 1:numel(tmp)));
tmpExpert = tmp(arrayfun(@(i) contains(tmp(i).name, '1_'), 1:numel(tmp)));
filesQA = arrayfun(@(f) strcat(DATADIR, f.name), tmpExpert, 'UniformOutput', false);
filesQANovice = arrayfun(@(f) strcat(DATADIR, f.name), tmpNovice, 'UniformOutput', false);


tmp = dir(DATADIR);
tmp = tmp(arrayfun(@(i) contains(tmp(i).name, 'out') & (~contains(tmp(i).name, 'QA')), 1:numel(tmp))); % remove stuff that is not the out files
tmpNovice = tmp(arrayfun(@(i) contains(tmp(i).name, '0_'), 1:numel(tmp)));
tmpExpert = tmp(arrayfun(@(i) contains(tmp(i).name, '1_'), 1:numel(tmp)));

areas = {'wS1', 'wS2', 'wM1', 'A1', 'V1', 'PPC', 'dCA1', 'mPFC', 'Striatum', 'wM2', 'ALM', 'tjM1'};

fields2keep = {'wPlain','ws','lambda','cvRsq','cvPart','yP','yTest','RsqST','MI','MIR','pRav','pWaldAv', ...
        'psWaldAv','PSTHgen','PSTHdata','PSTHRgen','cluster','n','mouse','date', ...
        'group','sess','sIndex','ClusterDepth','PeakToBaseline','MDS','ARAindex','struct','struct_acr','AP','ML'};
allResults.Neurons = [];


mergestruct = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);  % function used to merge structures

% Expert
if strcmp(VERSION, 'light')
for a=1:length(areas)
    disp(a)
    load(filesQA{a});
    area = fieldnames(QA);
    area = area{1};
    QAE.(area) = QA.(area);
    
    load(files{1});

    for i=1:length(results.Neurons)
        clear neuron
        for fi = 1:length(fields2keep)
            neuron.(fields2keep{fi}) = results.Neurons(i).(fields2keep{fi});
            neuron.area = results.area;
        end
        allResults.Neurons = cat(1, allResults.Neurons, neuron);
    end
end

GLM.allResultsE.Neurons = allResults.Neurons;
GLM.allResultsE.dm = results.Sessions(end).dm;

% Novice
allResults.Neurons = [];
for a=1:length(areas)
    disp(a)
    load(filesQANovice{a});
    area = fieldnames(QA);
    area = area{1};
    QAN.(area) = QA.(area);

    load(files{1})

    for i=1:length(results.Neurons)
        clear neuron
        for fi = 1:length(fields2keep)
            neuron.(fields2keep{fi}) = results.Neurons(i).(fields2keep{fi});
            neuron.area = results.area;
        end
        allResults.Neurons = cat(1, allResults.Neurons, neuron);
    end
end

GLM.allResultsN.Neurons = allResults.Neurons;
GLM.allResultsN.dm = results.Sessions(end).dm;
for a=1:length(areas)
    GLM.QA.(areas{a}) = mergestruct(QAE.(areas{a}), QAN.(areas{a}));
end

GLM.fitParams = fitParams;
GLM.fitOptions = fitOptions;
save(strcat(DATADIR, 'GLM.mat'), 'GLM', '-v7.3');
end  % end of light


if strcmp(VERSION, 'heavy')
% collect all results in 'GLMfulloutput' %
disp('Collecting the heavy structure ...')
for a=1:length(areas)
    disp(a)
    % experts
    load(filesQA{a});
    area = fieldnames(QA);
    area = area{1};
    QAE.(area) = QA.(area);
    
    load(files{1});

    GLMfulloutput.Expert.(area) = results;
    GLMfulloutput.Expert.fitOptions = fitOptions;
    GLMfulloutput.Expert.fitParams = fitParams;
    
    % novice
    load(filesQANovice{a});
    area = fieldnames(QA);
    area = area{1};
    QAN.(area) = QA.(area);

    load(files{1})

	
    GLMfulloutput.Novice.(area) = results;
    GLMfulloutput.Novice.fitOptions = fitOptions;
    GLMfulloutput.Novice.fitParams = fitParams;
end
save(strcat(DATADIR, 'GLMfulloutput.mat'), 'GLMfulloutput', '-v7.3');
end % end of heavy
