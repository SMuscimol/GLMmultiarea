%%% This file is meant to be called only from inside main.m %%%

%% collect all results in a lightweight structure 'GLM' %%
fields2keep = {'wPlain','ws','lambda','cvRsq','cvPart','yP','yTest','RsqST','MI','MIR','pRav','pWaldAv', ...
        'psWaldAv','PSTHgen','PSTHdata','PSTHRgen','cluster','n','mouse','date', ...
        'group','sess','sIndex','ClusterDepth','PeakToBaseline','MDS','ARAindex','struct','struct_acr','AP','ML'};
allResults.Neurons = [];

% Expert
for a=1:length(areas)
    disp(a)
    load(files{a})

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
    
    load(filesNovice{a})

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

GLM.QA = QA;

GLM.fitParams = fitParams;
GLM.fitOptions = fitOptions;
%% collect all results in 'GLMfulloutput' %%

for a=1:length(areas)
    area = areas{a};
    load(files{a})
    GLMfulloutput.Expert.(area) = results;
    GLMfulloutput.Expert.fitOptions = fitOptions;
    GLMfulloutput.Expert.fitParams = fitParams;
    load(filesNovice{a})
    GLMfulloutput.Novice.(area) = results;
    GLMfulloutput.Novice.fitOptions = fitOptions;
    GLMfulloutput.Novice.fitParams = fitParams;
end