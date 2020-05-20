function tIndexExpt = gettexptfromttest(tIndexTest, kFold, n, expt)

allIdxs.Hit = find(arrayfun(@(tidx) strcmp(expt.trial(tidx).tType,'Hit'),1:numel(expt.trial))); % indices of all hit trials in expt
allIdxs.CR = find(arrayfun(@(tidx) strcmp(expt.trial(tidx).tType,'CR'),1:numel(expt.trial))); % indices of all CR trials in expt


allExptIndices = sort([allIdxs.Hit(find(test(n.cvPart.Hit,kFold))), allIdxs.CR(find(test(n.cvPart.CR,kFold)))]); % all indices in expt used for this CV fold

tIndexExpt = allExptIndices(tIndexTest);

end